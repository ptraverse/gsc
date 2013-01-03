#! /usr/bin/env python
"""
jobs.py

Created by Lucas Swanson
Copyright (c) 2012 Canada's Michael Smith Genome Sciences Centre. All rights reserved.
"""

# import standard modules
import os

# import custom modules
from utils.log import GetLogPath
from utils.error import MyError
from utils.general import GetCommand
from utils.messages import LogMsg, DebugMsg
from utils.files_paths import (FileBoxCls, CheckFilePath,
  EnsureDirectoryExists)
from utils.subprocesses import (SortFile)

# handles all the jobs for all chromosomes
class HybridMetaJobCollectionCls: #{
  def __init__(self, options, log_info=None): #{
    self.options  = options
    self.log_info = log_info
    # one job collection per chromosome
    self.job_collections = dict()
  #} end def

  def __iter__(self): #{
    return self.job_collections.itervalues()
  #} end def

  def Chromosomes(self): #{
    return self.job_collections.iterkeys()
  #} end def

  def AddRecord(self, chromosome, record_str): #{
    if (chromosome not in self.job_collections): #{
      new_collection = HybridJobCollectionCls(chromosome, self.options)
      self.job_collections[chromosome] = new_collection
    #} end if
    self.job_collections[chromosome].AddRecord(record_str)
  #} end def
#} end class

# handles all the jobs for a single chromosome
class HybridJobCollectionCls: #{
  def __init__(self, chromosome, options): #{
    self.chromosome = chromosome
    # the important options are:
    #   "jobs_dir"
    #   "records_per_job"
    self.options = options
    self.CreateFullFile()
    self.jobs = list()
  #} end def

  def __iter__(self): #{
    return self.jobs.__iter__()
  #} end def

  def __del__(self): #{
    self.CloseFullFile()
  #} end def

  def CreateFullFile(self): #{
    EnsureDirectoryExists(self.options.jobs_dir)
    full_file_path = os.path.join(self.options.jobs_dir, "%s.%s.hyb.data" %
      (self.options.lib, self.chromosome))
    fail_msg = "cannot create chromosome %s all records file" % self.chromosome
    self.file = FileBoxCls(full_file_path, "w", fail_msg)
  #} end def

  def AddRecord(self, record_str): #{
    self.file.WriteLine(record_str)
  #} end def

  def SortFileByPosition(self): #{
    if (None == self.file): #{
      raise HybridJobError("all records file not created for "
        "chromosome %s" % self.chromosome)
    #} end if
    self.CloseFullFile()
    DebugMsg(self, "Sorting all records for chromosome %s" % self.chromosome)
    raise HybridJobError("not implemented! must choose sort column!")
    SortFile(self, self.file.path, 1, numeric=True)
  #} end def

  def SplitJobs(self): #{
    in_file = FileBoxCls(self.file.path, "r", "cannot read chromosome %s "
      "all records file" % self.chromosome)
    for record_str in in_file: #{
      if (0 == len(self.jobs) or self.jobs[-1].IsFull()): #{
        DebugMsg(self, "Adding new job to %s job collection, num jobs: %i" %
          (self.chromosome, len(self.jobs)+1))
        new_job = HybridJobCls(self.chromosome, len(self.jobs)+1,
          self.options, log_info=self.log_info)
        self.jobs.append(new_job)
      #} end if
      self.jobs[-1].AddRecord(record_str)
    #} end for
  #} end def

  def CloseFullFile(self): #{
    if (hasattr(self, "file") and None != self.file): #{
      self.file.Close()
    #} end if
  #} end def
#} end class

# a single job, all records for the same chromosome,
#   and sorted by chromosome position
class HybridJobCls: #{
  def __init__(self, chromosome, index, options, log_info=None): #{
    self.chromosome  = chromosome
    self.index       = index
    self.num_records = 0
    self.file        = None
    self.options     = options
    self.log_info    = log_info
  #} end def

  def IsFull(self): #{
    if (self.num_records < self.options.records_per_job): #{
      return False
    #} end if
    return True
  #} end def

  def CreateJobDataFile(self): #{
    job_dir = os.path.join(self.options.jobs_dir,
      "job_%i" % self.index)
    EnsureDirectoryExists(job_dir)
    job_path = os.path.join(job_dir, "%s.%s.%i.hyb.data" %
      (self.options.lib, self.chromosome, self.index))
    fail_msg = "cannot create job data file"
    self.file = FileBoxCls(job_path, "w", fail_msg)
  #} end def

  def AddRecord(self, record_str): #{
    if (self.IsFull()): #{
      raise HybridJobError("cannot add record to job %s.%i, job is "
        "full\n  %s" % (self.chromosome, self.index, record_str))
    #} end if
    if (None == self.file): #{
      self.CreateJobDataFile()
    #} end if
    self.file.WriteLine(record_str)
  #} end def

  def Command(self, local=False): #{
    python_path = GetCommand(self, "python", local)
    script_path = os.path.join(os.environ['BARNACLE_PATH'],
      "support", "hybrid", "calculate.py")
    CheckFilePath(script_path, "hybrid read-support script")
    command_list = ["time", python_path, script_path, self.options.lib,
      self.file.path, self.options.p2g_path, self.options.r2c_path,
      "--min-overlap", "%i" % self.options.min_overlap]
    if (self.log_info['debug']): #{
      command_list.append("--debug")
    #} end if
    command = " ".join(command_list)
    if (local): #{
      if ("hyb" in os.path.basename(self.file.path)): #{
        log_type = "calculate"
      else:
        log_type = "hyb.calculate"
      #} end if
      log_path = GetLogPath(self.file.path, log_type, options.output_dir)
      command += " --log-file %s" % log_path
      return command
    #} end if
    full_command_list = [
      "export BARNACLE_PATH=%s" % os.environ["BARNACLE_PATH"],
      "export PATH=$BARNACLE_PATH:$PATH",
      "export PYTHONPATH=.:$BARNACLE_PATH:$PYTHONPATH",
      command,
    ]
    full_command = ";".join(full_command_list)
    return full_command
  #} end def

  def CommandPath(self): #{
    return os.path.join(self.options.jobs_dir,
      "%s-hyb.%s.%i.sh" % (self.options.lib, self.chromosome, self.index))
  #} end def

  def StdOutPath(self): #{
    return "%s.o1" % self.CommandPath()
  #} end def

  def StdErrPath(self): #{
    return "%s.e1" % self.CommandPath()
  #} end def

  def Close(self): #{
    if (None != self.file): #{
      self.file.Close()
    #} end if
  #} end def
#} end class

#### EXCEPTION CLASSES ####
class HybridJobError(MyError): #{
  pass
#} end class
