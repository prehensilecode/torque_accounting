#!/usr/bin/env python

# Author: David Chin <dwchin . acm.org>
# License: GPLv3

import sys, os, re
import datetime
import copy
import cPickle
import numpy as np

from optparse import OptionParser

# States are: 
# - Q - queued
# - S - started
# - D - deleted
# - E - ended
#
# Possible state transitions:
# - normal termination: Q -> S -> E
# - job qdel'ed: Q -> S -> D -> E
#
# Each log line is a possible new state for a job, with a timestamp
#
# Each state has possibly different data in the line:
# 
# Q: 01/01/2013 20:40:57;Q;134907.rhel6pbs.deac.wfu.edu;queue=rhel6
# - this is the least useful state since it only gives jobid as, e.g. 123456[].server.example.com
#
# S: 01/01/2013 00:18:53;S;134394.rhel6pbs.deac.wfu.edu;user=ssajuthi group=langefeldGrp jobname=impute.chr18.15.pbs queue=rhel6 ctime=1357014926 qtime=1357014926 etime=1357014926 start=1357017533 owner=ssajuthi@rhel6head1.deac.wfu.edu exec_host=bcg01bl11/7 Resource_List.arch=x86_64 Resource_List.cput=96:00:00 Resource_List.mem=10000mb Resource_List.ncpus=1 Resource_List.neednodes=1:ppn=1:ethernet Resource_List.nodect=1 Resource_List.nodes=1:ppn=1:ethernet Resource_List.pmem=5000mb Resource_List.walltime=96:00:00 
# 
# D: 01/01/2013 17:20:10;D;134596.rhel6pbs.deac.wfu.edu;requestor=leplnd6@rhel6head1.deac.wfu.edu
#
# E: 01/01/2013 17:20:11;E;134596.rhel6pbs.deac.wfu.edu;user=leplnd6 group=physp2Grp jobname=Li3PO4beta queue=rhel6 ctime=1357078656 qtime=1357078656 etime=1357078656 start=1357078657 owner=leplnd6@rhel6head1.deac.wfu.edu exec_host=bc07bl03/7+bc07bl03/6+bc07bl03/5+bc07bl03/4+bc07bl03/3+bc07bl03/2+bc07bl03/1+bc07bl03/0+bc07bl04/7+bc07bl04/6+bc07bl04/5+bc07bl04/4+bc07bl04/3+bc07bl04/2+bc07bl04/1+bc07bl04/0 Resource_List.arch=x86_64 Resource_List.cput=800:00:00 Resource_List.mem=45gb Resource_List.neednodes=2:ppn=8:ethernet Resource_List.nodect=2 Resource_List.nodes=2:ppn=8:ethernet Resource_List.pmem=1024mb Resource_List.walltime=100:00:00 session=10413 end=1357078811 Exit_status=271 resources_used.cput=00:40:00 resources_used.mem=24929396kb resources_used.vmem=32561732kb resources_used.walltime=00:02:34

class Memory:
    """
    Represents amount of memory
    """
    def __init__(self, mem=None):
        # Example inputs:
        # - '1024mb'
        # - {'qty': 4239.123, 'units': 'GiB'}
        # - 123
        # - 123.456
        self.__KILO = 1024.
        self.__MEGA = 1048576.
        self.__kiB = 'kiB'
        self.__MiB = 'MiB'
        self.__GiB = 'GiB'

        if mem:
            if type(mem) == str:
                self.__mem = self.__convert_memory(mem, 'kb')
            elif type(mem) == dict:
                self.__mem = mem
            elif type(mem) == int or type(mem) == float:
                # assume kiB
                memstr = ''.join([str(mem), 'kb'])
                self.__mem = {'qty': mem, 'units': self.__kiB}
            else:
                raise Exception
        else:
            self.__mem = {'qty': 0., 'units': self.__kiB}

        self.qty = self.__mem['qty']
        self.units = self.__mem['units']

        self.__to_kiB()


    def in_MiB(self):
        """Returns float of qty in MiB"""
        qty = 0.
        if self.units == self.__kiB:
            qty = self.qty / self.__KILO
        elif self.units == self.__GiB:
            qty = self.qty * self.__KILO

        return qty


    def in_GiB(self):
        """Returns float of qty in GiB"""
        qty = 0.
        if self.units == self.__kiB:
            qty = self.qty / self.__MEGA
        elif self.units == self.__MiB:
            qty = self.qty / self.__KILO

        return qty


    def str_in_MiB(self):
        qty = self.in_MiB()
        units = self.__MiB
        formatstr = ""
        if qty < 1:
            formatstr = "{qty:.2e} {units:3.3}"
        else:
            formatstr = "{qty:.2f} {units:3.3}"

        return formatstr.format(qty=qty, units=units)


    def str_in_GiB(self):
        qty = self.in_GiB()
        units = self.__GiB
        if qty < 1:
            formatstr = "{qty:.2e} {units:3.3}"
        else:
            formatstr = "{qty:.2f} {units:3.3}"
        return formatstr.format(qty=qty, units=units)


    def copy(self):
        return copy.deepcopy(self)


    def pretty_print(self):
        """Select appropriate units based on quantity"""
        qty = 0.
        units = ''
        if self.__mem['qty'] < 1:
            qty = self.__mem['qty'] * self.__KILO
            units = 'B'
        elif self.in_MiB() < 1:
            qty = self.__mem['qty']
            units = self.__mem['units']
        elif self.in_GiB() < 1:
            qty = self.in_MiB()
            units = self.__MiB
        else:
            qty = self.in_GiB()
            units = self.__GiB

        formatstr = ""
        if qty < 1:
            formatstr = "{qty:.2e} {units:<}"
        else:
            formatstr = "{qty:.2f} {units:<}"
        return formatstr.format(qty=qty, units=units)
        

    def __to_kiB(self):
        """convert to kiB"""

        if self.__mem['units'] != self.__kiB:
            if self.__mem['units'] == self.__MiB:
                self.__mem['qty'] *= self.__KILO
            elif self.__mem['units'] == self.__GiB:
                self.__mem['qty'] *= self.__MEGA

            self.__mem['units'] = self.__kiB

            self.qty = self.__mem['qty']
            self.units = self.__mem['units']


    def __convert_memory(self, memstr=None, units=None):
        """
        Converts TORQUE memory strings, in the form, 'NNNNNNkb' to dict: {'qty': QTY, 'units': UNITS}
        """
        memkb = 0.

        bpat = re.compile(r'(\d+\.?\d*)(b)$', re.I)
        kbpat = re.compile(r'(\d+\.?\d*)(ki?b)$', re.I)
        mbpat = re.compile(r'(\d+\.?\d*)(mi?b)$', re.I)
        gbpat = re.compile(r'(\d+\.?\d*)(gi?b)$', re.I)

        if memstr:
            bs = bpat.search(memstr)
            ks = kbpat.search(memstr)
            ms = mbpat.search(memstr)
            gs = gbpat.search(memstr)

            if bs:
                memkb = float(memstr.split(bs.group(2))[0]) / self.__KILO
            elif kbpat.match(memstr):
                memkb = float(memstr.split(ks.group(2))[0])
            elif mbpat.match(memstr):
                memkb = float(memstr.split(ms.group(2))[0]) * self.__KILO
            elif gbpat.match(memstr):
                memkb = float(memstr.split(gs.group(2))[0]) * self.__MEGA
            else:
                sys.stderr.write('convert_memory(): cannot parse "{0}"\n'.format(memstr))
                raise Exception

        if not units:
            if memkb < self.__KILO:
                units = 'kB'
            elif memkb < self.__MEGA:
                units = 'MB'
            else:
                units = 'GB'

        mem = {}
        if units == 'GB' or units == 'GiB' or units == 'gb':
            mem['units'] = 'GiB'
            mem['qty'] = memkb/(self.__MEGA)
        elif units == 'MB' or units == 'MiB' or units == 'mb':
            mem['units'] = 'MiB'
            mem['qty'] = memkb/(self.__KILO)
        elif units == 'kB' or units == 'kiB' or units == 'kb':
            mem['units'] = 'kiB'
            mem['qty'] = memkb
        return mem


    def __cmp__(self, other):
        return (self.__mem['qty'] - other.__mem['qty'])

    
    def __add__(self, other):
        return Memory(self.__mem['qty'] + other.__mem['qty'])


    def __sub__(self, other):
        return Memory(self.__mem['qty'] - other.__mem['qty'])

    
    def __repr__(self):
        formatstr = "{qty:.2f} {units:3.3}"
        if self.qty < 1:
            formatstr = "{qty:.2e} {units:3.3}"
        else:
            formatstr = "{qty:.2f} {units:3.3}"
        return formatstr.format(qty=self.qty, units=self.units)


    def __str__(self):
        return self.__repr__()
        

class JobEvent:
    def __init__(self, event=None, timestamp=None):
        self.event = event
        self.timestamp = self.__parse_timestamp(timestamp)

    def __parse_timestamp(self, timestamp):
        thisdate, thistime = timestamp.split(' ')
        date = [int(d) for d in thisdate.split('/')]
        timeofday = [int(t) for t in thistime.split(':')]
        return (datetime.datetime(date[2], date[0], date[1], timeofday[0], timeofday[1], timeofday[2]))

    def copy(self):
        return copy.deepcopy(self)

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        formatstr = "Timestamp: {timestamp}    Event: {event}"
        return formatstr.format(timestamp=self.timestamp, event=self.event)


class Job:
    def __init__(self, jobid=None, event=None, timestamp=None, joblog=None):
        self.__debug_p = False

        self.jobid = jobid
        self.event_list = []
        self.event_list.append(JobEvent(event, timestamp))

        self.queue = None
        self.user = None
        self.group = None
        self.jobname = None
        self.submithost = None
        self.qdel_requestor = None

        self.owner = None

        self.ctime = None
        self.qtime = None
        self.etime = None

        self.start = None
        self.end = None
        self.exit_status = None

        self.exec_host = None

        self.ncpus = None
        self.unique_hosts = None
        self.nnodes = None

        self.resources_requested = {}
        self.resources_used = {}

        self.__parse_joblog(joblog)


    def new_event(self, event=None, timestamp=None, joblog=None):
        """Insert a new event record for the job"""
        self.event_list.append(JobEvent(event, timestamp))

        self.__parse_joblog(joblog)


    def get_latest_event(self):
        """Latest event recorded for the job"""
        return self.event_list[-1]


    def duration(self):
        """Amount of walltime the job ran"""
        retval = None
        if 'walltime' in self.resources_used:
            retval = self.resources_used['walltime']
        elif self.start and self.end:
            retval = self.end - self.start
        return retval


    def wait_time(self):
        """Amount of time the job waited in the queue before executing"""
        retval = None
        if self.start and self.end:
            retval = self.start - self.qtime
        return retval
            

    def walltime_discrepancy(self):
        """(Req. walltime - Actual walltime)/Actual walltime"""
        retval = None
        d = self.duration()
        if d:
            retval = (self.resources_requested['walltime'].total_seconds() - d.total_seconds())/d.total_seconds()
        return retval

    def is_complete(self):
        """If job has a complete record, i.e. has a Start event and an End event"""
        has_s = False
        has_e = False

        for ev in self.event_list:
            if ev.event == 'S':
                has_s = True
            elif ev.event == 'E':
                has_e = True

        return (has_s and has_e)

    
    def __convert_timedelta_str(self, timedelta_str):
        td = [int(t) for t in timedelta_str.split(':')]
        return datetime.timedelta(hours=td[0], minutes=td[1], seconds=td[2])


    def __set_member_maybe(self, membername, data):
        if self.__debug_p:
            print 'DEBUG: membername = ', membername
            print 'DEBUG: data[', membername, '] = ', data[membername]

        if membername in data:
            reslist_pat = re.compile(r'^Resource_List\.(\w+)', re.I)
            resused_pat = re.compile(r'^resources_used\.(\w+)', re.I)

            reslist_search = reslist_pat.search(membername)
            resused_search = resused_pat.search(membername)

            if reslist_search:
                #resname = reslist_search.group(1)
                resname = membername.split('.')[1]
                if resname == 'mem' or resname == 'pmem':
                    self.resources_requested[resname] = Memory(data[membername])
                elif resname == 'cput' or resname == 'walltime':
                    self.resources_requested[resname] = self.__convert_timedelta_str(data[membername])
                elif resname == 'nodect':
                    self.resources_requested[resname] = int(data[membername])
                else:
                    self.resources_requested[resname] = data[membername]
            elif resused_search:
                #resname = resused_search.group(1)
                resname = membername.split('.')[1]
                if resname == 'mem' or resname == 'vmem':
                    self.resources_used[resname] = Memory(data[membername])
                elif resname == 'cput' or resname == 'walltime':
                    self.resources_used[resname] = self.__convert_timedelta_str(data[membername])
                else:
                    self.resources_used[resname] = data[membername]
            elif membername == 'owner':
                # technically not needed since owner is always user@submithost
                if not self.__dict__[membername]:
                    self.owner = data[membername].split('@')[0]
                if not self.__dict__['submithost']:
                    self.submithost = data[membername].split('@')[1]
            elif membername == 'ctime' or membername == 'qtime' or membername == 'etime' or membername == 'start' or membername == 'end':
                if not self.__dict__[membername]:
                    self.__dict__[membername] = datetime.datetime.fromtimestamp(int(data[membername]))
            elif membername == 'exec_host':
                if not self.__dict__[membername]:
                    self.__parse_exec_host(data[membername])
            elif membername == 'Exit_status':
                if not self.__dict__[membername.lower()]:
                    self.__dict__[membername.lower()] = int(data[membername])
            elif not self.__dict__[membername.lower()]:
                self.__dict__[membername.lower()] = data[membername]


    def __parse_exec_host(self, exec_host_str):
        node_cpu = exec_host_str.split('+')
        self.ncpus = len(node_cpu)

        self.unique_hosts = set([])
        for nc in node_cpu:
            node = nc.split('/')[0]
            self.unique_hosts.add(node)

        self.nnodes = len(self.unique_hosts)


    def __parse_joblog(self, joblog):
        if self.__debug_p:
            print 'DEBUG: joblog=', joblog
        cur_event = self.get_latest_event().event
        joblogdata = joblog.split(' ')
        data = {}
        for jd in joblogdata:
            dataelem = jd.split('=')
            data[dataelem[0]] = dataelem[1]

        ### The Q state gives no useful information except for the queue
        ### and the timestamp. But the timestamp is the same as qtime
        if cur_event == 'Q':
            # joblog gives only the queue name: queue=rhel6
            self.queue = data['queue']
        elif cur_event == 'S':
            self.__set_member_maybe('queue', data)
            self.__set_member_maybe('user', data)
            self.__set_member_maybe('group', data)
            self.__set_member_maybe('jobname', data)
            self.__set_member_maybe('owner', data)
            self.__set_member_maybe('ctime', data)
            self.__set_member_maybe('qtime', data)
            self.__set_member_maybe('etime', data)
            self.__set_member_maybe('start', data)
            self.__set_member_maybe('exec_host', data)
            self.__set_member_maybe('Resource_List.mem', data)
            self.__set_member_maybe('Resource_List.cput', data)
            self.__set_member_maybe('Resource_List.walltime', data)
            self.__set_member_maybe('Resource_List.nodect', data)
        elif cur_event == 'D':
            self.qdel_requestor = data['requestor']
        elif cur_event == 'E':
            self.__set_member_maybe('queue', data)
            self.__set_member_maybe('user', data)
            self.__set_member_maybe('group', data)
            self.__set_member_maybe('jobname', data)
            self.__set_member_maybe('owner', data)
            self.__set_member_maybe('ctime', data)
            self.__set_member_maybe('qtime', data)
            self.__set_member_maybe('etime', data)
            self.__set_member_maybe('start', data)
            self.__set_member_maybe('end', data)
            self.__set_member_maybe('exec_host', data)
            self.__set_member_maybe('Exit_status', data)
            self.__set_member_maybe('Resource_List.mem', data)
            self.__set_member_maybe('Resource_List.cput', data)
            self.__set_member_maybe('Resource_List.walltime', data)
            self.__set_member_maybe('Resource_List.nodect', data)
            self.__set_member_maybe('resources_used.mem', data)
            self.__set_member_maybe('resources_used.vmem', data)
            self.__set_member_maybe('resources_used.cput', data)
            self.__set_member_maybe('resources_used.walltime', data)
        else:
            print 'UNKNOWN EVENT:', cur_event
        
    def printout(self):
        for k,v in sorted(self.__dict__.iteritems()):
            print "{k}: {v}".format(k=k, v=v)
        print "duration:", self.duration()

    def __str__(self):
        if self.owner != self.user:
            print 'FOOBAR! owner != user'
        retstr = ''
        if self.qdel_requestor:
            formatstr = "Job ID: {self.jobid}\n    Jobname: {self.jobname}\n    Queue: {self.queue}\n    User: {self.user}\n    Group: {self.group}\n    Submit host: {self.submithost}\n    Event: {event}\n    Qdel requestor: {self.qdel_requestor}\n    Is complete: {comp}"
            retstr = formatstr.format(self=self, event=[str(s) for s in self.event_list], comp=self.is_complete())
        else:
            formatstr = "Job ID: {self.jobid}\n    Jobname: {self.jobname}\n    Queue: {self.queue}\n    User: {self.user}\n    Group: {self.group}\n    Submit host: {self.submithost}\n    Event: {event}\n    Is complete: {comp}"
            retstr = formatstr.format(self=self, event=[str(s) for s in self.event_list], comp=self.is_complete())
        return retstr



def group_stats(job_dict, groupname):
    thres = datetime.timedelta(seconds=59)
    n_jobs = 0
    wait_times = []
    durations = []
    for jobid,job in sorted(job_dict.iteritems()):
        if job.group == groupname:
            #print jobid, job.user, job.group, job.wait_time(), job.duration()
            wt = job.wait_time()
            du = job.duration()
            if wt and du:
                if du > thres:
                    wait_times.append(wt.total_seconds())
                    durations.append(du.total_seconds())
                    n_jobs += 1

    wt_arr = np.array(wait_times, np.uint64)
    du_arr = np.array(durations, np.uint64)

    print "Stats for group {g}".format(g=groupname)
    print "{n} jobs of duration > {th} found".format(n=n_jobs, th=thres)
    print ""

    if n_jobs:
        print "WAITING TIMES:"
        print "Min: ", datetime.timedelta(seconds=float(wt_arr.min()))
        print "Max: ", datetime.timedelta(seconds=float(wt_arr.max()))
        print "Mean:", datetime.timedelta(seconds=float(wt_arr.mean()))

        print ""
        print "DURATIONS:"
        print "Min: ", datetime.timedelta(seconds=float(du_arr.min()))
        print "Max: ", datetime.timedelta(seconds=float(du_arr.max()))
        print "Mean:", datetime.timedelta(seconds=float(du_arr.mean()))


def user_stats(job_dict, username):
    thres = datetime.timedelta(seconds=59)
    n_jobs = 0
    wait_times = []
    durations = []
    for jobid,job in sorted(job_dict.iteritems()):
        if job.user == username:
            wt = job.wait_time()
            du = job.duration()
            if wt and du:
                if du > thres:
                    wait_times.append(wt.total_seconds())
                    durations.append(du.total_seconds())
                    n_jobs += 1

    wt_arr = np.array(wait_times, np.uint64)
    du_arr = np.array(durations, np.uint64)

    print "Stats for user {u}".format(u=username)
    print "{n} jobs of duration > {th} found".format(n=n_jobs, th=thres)
    print ""

    if n_jobs:
        print "WAITING TIMES:"
        print "Min: ", datetime.timedelta(seconds=float(wt_arr.min()))
        print "Max: ", datetime.timedelta(seconds=float(wt_arr.max()))
        print "Mean:", datetime.timedelta(seconds=float(wt_arr.mean()))

        print ""
        print "DURATIONS:"
        print "Min: ", datetime.timedelta(seconds=float(du_arr.min()))
        print "Max: ", datetime.timedelta(seconds=float(du_arr.max()))
        print "Mean:", datetime.timedelta(seconds=float(du_arr.mean()))


def walltime_stats(job_dict):
    n_jobs = 0
    durations = []
    discrepancies = []
    thres = datetime.timedelta(seconds=59)
    for jobid,job in sorted(job_dict.iteritems()):
        d = job.duration()
        if d and d > thres:
            n_jobs += 1
            durations.append(d.total_seconds())
            discrepancies.append(job.walltime_discrepancy())

    du_arr = np.array(durations, np.uint64)
    ds_arr = np.array(discrepancies, np.float64)

    print "WALLTIME STATS (actual)"
    print "{n} jobs found (duration > {t})".format(n=n_jobs, t=thres)
    print "discrepancy = (req. walltime - actual walltime)/actual walltime"
    print ""

    if n_jobs:
        print "Min. walltime: ", datetime.timedelta(seconds=float(du_arr.min()))
        print "Max. walltime: ", datetime.timedelta(seconds=float(du_arr.max()))
        print "Mean walltime: ", datetime.timedelta(seconds=float(du_arr.mean()))

        print ""

        print "Min. discrepancy:", datetime.timedelta(seconds=float(ds_arr.min()))
        print "Max. discrepancy:", datetime.timedelta(seconds=float(ds_arr.max()))
        print "Mean discrepancy:", datetime.timedelta(seconds=float(ds_arr.mean()))


def parse_line(line):
    """Returns (timestamp, event, jobid, accounting record)"""
    thisline = line.strip().split(';')
    return (thisline[0], thisline[1], thisline[2], thisline[3])

def main(opt, args):
    job_dict = {}
    for a in args:
        try:
            f = open(a, 'ro')
            try:
                lines = f.readlines()
            except IOError:
                print "ERROR:"
        finally:
            f.close()

        #print("Got {n} lines".format(n=len(lines)))

        nlines = len(lines)
        for i in range(nlines):
            timestamp, event, jobid, joblog = parse_line(lines[i])
            if not jobid in job_dict:
                job_dict[jobid] = Job(jobid=jobid, event=event, timestamp=timestamp, joblog=joblog)
            else:
                job_dict[jobid].new_event(event=event, timestamp=timestamp, joblog=joblog)



    # Fix up job array jobs
    # Job array jobs have only one 'Q' entry, which bears a 
    # generic job ID like 123456[].rhel6pbs.deac.wfu.edu
    # Want to add 'Q' events to all the actual array jobs, 
    # and then delete that generic initial Q event
    jobarray_base_ids = []
    jobarraybasepat = re.compile(r'^(\d+)\[\]\.\w*')
    for k,v in job_dict.iteritems():
        jabs = jobarraybasepat.search(k)
        if jabs:
            jobarray_base_ids.append(jabs.group(1))

    jobarraypat = re.compile(r'^(\d+)\[\d+\]\.[\w\.]+')
    for k,v in job_dict.iteritems():
        jas = jobarraypat.search(k)
        if jas:
            for jid in jobarray_base_ids:
                if jid == jas.group(1):
                    basejobid = ''.join([jid, '[].rhel6pbs.deac.wfu.edu'])
                    base_event = job_dict[basejobid].event_list[0].copy()
                    v.event_list.insert(0, base_event)

    for i in jobarray_base_ids:
        del job_dict[''.join([i, '[].rhel6pbs.deac.wfu.edu'])]
    
    print("Found {n} jobs".format(n=len(job_dict)))

    ncomplete = 0
    has_duration = 0
    td0 = datetime.timedelta(seconds=59)
    for k,v in job_dict.iteritems():
        if v.is_complete():
            ncomplete += 1
        if v.duration() and v.duration() > td0:
            has_duration += 1
        v.printout()
        if v.resources_used:
            for res in ['walltime', 'cput', 'vmem']:
                print "resources_used[{res}]: {val}".format(res=res, val=v.resources_used[res])
            print "walltime discrepancy:", v.walltime_discrepancy()
        print ""

    print("Found {n} jobs with complete log record".format(n=ncomplete))
    print("Found {n} with duration > {td0}".format(n=has_duration, td0=td0))


    #durations = []
    #for k,v in job_dict.iteritems():
    #    if v.duration() and v.duration() > td0:
    #        durations.append(v.duration().total_seconds() * 1000000)

    #dur = np.array(durations, dtype=np.uint64)
    
    #print 'Max:', datetime.timedelta(seconds=(dur.max()/1000000.))
    #print 'Min:', datetime.timedelta(seconds=(dur.min()/1000000.))
    #print 'Mean:', datetime.timedelta(seconds=(dur.mean()/1000000.))

    #print ''
    #print ''
    print ''

    group_stats(job_dict, 'langefeldGrp')

    print("-----------------------------------------------------------------")

    group_stats(job_dict, 'thonhauserGrp')

    print("-----------------------------------------------------------------")

    group_stats(job_dict, 'salsburyGrp')

    print("-----------------------------------------------------------------")

    user_stats(job_dict, 'canepap')

    print("-----------------------------------------------------------------")

    user_stats(job_dict, 'negureal')

    print("-----------------------------------------------------------------")

    walltime_stats(job_dict)


if __name__ == '__main__':
    parser = OptionParser()
    (opt, args) = parser.parse_args()

    main(opt, args)

