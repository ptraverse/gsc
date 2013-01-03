#! /usr/bin/env python
"""
subprocesses.py

Created by Lucas Swanson
Copyright (c) 2012 Canada's Michael Smith Genome Sciences Centre. All rights reserved.
"""

# import standard modules
import os, shlex, signal, subprocess

# import custom modules
from error import MyError
from messages import ErrMsg, DebugMsg
from files_paths import OpenFile

# CONSTANTS
STRING_IN  = subprocess.PIPE              # -1
STRING_OUT = subprocess.PIPE              # -1
STD_OUT    = subprocess.STDOUT            # -2
STREAM_OUT = min(STRING_OUT, STD_OUT) - 1 # -3

def OpenCommandFile(file_info, default_mode): #{
  if (file_info in [STRING_OUT, STD_OUT, STREAM_OUT] or
      isinstance(file_info, file)): #{
    return file_info
  elif (None != file_info):
    # if just a path was passed in, use default values
    if (isinstance(file_info, str)): #{
      file_info = {'path': file_info, 'mode':default_mode, 'fail_msg':None}
    #} end if
    return OpenFile \
      (file_info['path'], file_info['mode'], file_info['fail_msg'])
  #} end if
  return None
#} end def

def OpenCommandFiles(stdin, stdout, stderr): #{
  stdin_file = OpenCommandFile(stdin, "r")
  stdout_file = OpenCommandFile(stdout, "w")
  stderr_file = OpenCommandFile(stderr, "w")
  return (stdin_file, stdout_file, stderr_file)
#} end def

def CloseCommandFile(file): #{
  if (None != file and file not in [STRING_OUT, STD_OUT, STREAM_OUT]): #{
    file.close()
  #} end if
# end def

def CloseCommandFiles(stdin_file, stdout_file, stderr_file): #{
  CloseCommandFile(stdin_file)
  CloseCommandFile(stdout_file)
  CloseCommandFile(stderr_file)
#} end def

# dpt = disable profiling timer flag
def RunCommandFromList(cmd_list, stdin=None, stdout=None, stderr=None,
    input=None, dpt=False): #{
  if (None != input and (STRING_IN != stdin or STRING_OUT != stdout)): #{
    ErrMsg("Warning RunCommandFromList only uses input if stdin is "
        "STRING_IN and stdout is STRING_OUT.")
  #} end if
  stdout_result = stderr_result = result = None
  # open files
  (stdin_file, stdout_file, stderr_file) = OpenCommandFiles(
    stdin, stdout, stderr)
  try:
    # disable profiling timer
    if (dpt): #{
      try:
        (old_delay, old_interval) = signal.setitimer(signal.ITIMER_PROF, 0)
      except AttributeError, e:
        ErrMsg("Python version does not support disabling profiling timer")
        dpt = False
      # end try
    #} end if
    try:
      # run command
      if (STREAM_OUT == stdout_file): #{
        #ErrMsg("Using Popen().stdout")
        proc = subprocess.Popen(cmd_list, shell=False, cwd=os.getcwd(),
          stdin=stdin_file, stdout=STRING_OUT, stderr=stderr_file)
        stdout_result = proc.stdout
      elif (STRING_OUT == stdout_file or STRING_OUT == stderr_file):
        #ErrMsg("Using Popen().communicate()")
        # use Popen() and communicate() to save the stdout and stderr output
        proc = subprocess.Popen(cmd_list, shell=False, cwd=os.getcwd(),
          stdin=stdin_file, stdout=stdout_file, stderr=stderr_file)
        (stdout_result, stderr_result) = proc.communicate(input)
        #ErrMsg("stdout: %s\nstderr: %s" % (stdout_result, stderr_result))
        if (0 != proc.returncode): #{
          raise MyError("error running command: %s\n%i" %
            (" ".join(cmd_list), proc.returncode))
        #} end if
      else:
        # just use call() to run the command
        #ErrMsg("Using call()")
        result = subprocess.call(cmd_list, shell=False, cwd=os.getcwd(),
          stdin=stdin_file, stdout=stdout_file, stderr=stderr_file)
      #} end if
    except OSError, e:
      raise MyError("error running command: %s\n%s" %
        (" ".join(cmd_list), e))
    finally:
      # re-enable profiling timer
      if (dpt): #{
        signal.setitimer(signal.ITIMER_PROF, old_delay, old_interval)
      #} end if
    # end try
  finally:
    # close files
    CloseCommandFiles(stdin_file, stdout_file, stderr_file)
  # end try
  if (STREAM_OUT == stdout_file): #{
    #ErrMsg("Returning stdout result as stream")
    return stdout_result
  #} end if
  if (STRING_OUT == stdout_file or STRING_OUT == stderr_file): #{
    #ErrMsg("Returning stdout and stderr results")
    #ErrMsg("stdout: %s" % stdout_result)
    #ErrMsg("stderr: %s" % stderr_result)
    return (stdout_result, stderr_result)
  #} end if
  #ErrMsg("Returning call() result")
  return result
#} end def

# dpt = disable profiling timer flag
def RunCommandFromString(cmd_string, stdin=None, stdout=None, stderr=None,
    input=None, dpt=False):
  cmd_list = shlex.split(cmd_string)
  return RunCommandFromList(cmd_list, stdin=stdin, stdout=stdout,
    stderr=stderr, input=input, dpt=dpt)
#} end def

# dpt = disable profiling timer flag
def RunPipeFromLists(cmd_lists, stdin=None, stdout=None, stderr=None,
    dpt=False):
  if (1 == len(cmd_lists)): #{
    return RunCommandFromList(cmd_lists[0], stdin=stdin, stdout=stdout,
      stderr=stderr, dpt=dpt)
  #} end if
  stdout_result = stderr_result = result = None
  # open files
  (stdin_file, stdout_file, stderr_file) = \
    OpenCommandFiles(stdin, stdout, stderr)
  try:
    # disable profiling timer
    if (dpt): #{
      try:
        (old_delay, old_interval) = signal.setitimer(signal.ITIMER_PROF, 0)
      except AttributeError, e:
        ErrMsg("Python version does not support disabling profiling timer")
        dpt = False
      # end try
    #} end if
    try:
      processes = list()
      # open processes for pipe
      for (process_index, cmd_list) in enumerate(cmd_lists): #{
        # if this is the first process, use the stdin_file
        if (0 == process_index): #{
          curr_stdin_file = stdin_file
        # else use the previous process's output
        else:
          curr_stdin_file = processes[process_index - 1].stdout
          # Allow p1 to receive a SIGPIPE if p2 exits.
          processes[process_index - 1].stdout.close()
        #} end if
        processes.append(subprocess.Popen(cmd_list,
          stdin=curr_stdin_file, stdout=subprocess.PIPE,
          preexec_fn=lambda: signal.signal(signal.SIGPIPE, signal.SIG_DFL)))
      #} end for
      # run final command
      if (STRING_OUT == stdout_file or STRING_OUT == stderr_file): #{
        #ErrMsg("Using Popen()")
        # use Popen() and communicate() to save the stdout and stderr output
        (stdout_result, stderr_result) = subprocess.Popen(cmd_list,
          stdin=processes[-1].stdout, stdout=stdout_file,
          stderr=stderr_file).communicate()
        #ErrMsg("stdout: %s\nstderr: %s" % (stdout_result, stderr_result))
      else:
        # just use call() to run the command
        #ErrMsg("Using call()")
        result = subprocess.call(cmd_list,
          stdin=processes[-1].stdout, stdout=stdout_file, stderr=stderr_file)
      #} end if
    finally:
      # re-enable profiling timer
      if (dpt): #{
        signal.setitimer(signal.ITIMER_PROF, old_delay, old_interval)
      #} end if
    # end try
  finally:
    # close files
    CloseCommandFiles(stdin_file, stdout_file, stderr_file)
  # end try
  if (STRING_OUT == stdout_file or STRING_OUT == stderr_file): #{
    #ErrMsg("Returning results")
    #ErrMsg("stdout: %s" % stdout_result)
    #ErrMsg("stderr: %s" % stderr_result)
    return (stdout_result, stderr_result)
  #} end if
  #ErrMsg("Returning call() result")
  return result
#} end def

# dpt = disable profiling timer flag
def RunPipeFromStrings(cmd_strings, stdin=None, stdout=None, stderr=None,
    dpt=False):
  if (1 == len(cmd_strings)): #{
    return RunCommandFromString(cmd_strings[0], stdin=stdin, stdout=stdout,
      stderr=stderr, dpt=dpt)
  #} end if
  cmd_lists = [shlex.split(cmd_string) for cmd_string in cmd_strings]
  return RunPipeFromLists(cmd_lists, stdin=stdin, stdout=stdout,
    stderr=stderr, dpt=dpt)
#} end def

def SortFile(self, path, sort_field, numeric=False): #{
  if (numeric): #{
    sort_field = "%in" % sort_field
  else:
    sort_field = "%i" % sort_field
    # set the LC_ALL environment variable to ensure sort uses byte values
    os.environ["LC_ALL"] = "C"
  #} end if
  sort_cmd_list = [ "sort", "--key", sort_field, "--temporary-directory",
    os.path.dirname(path), "--output", path, path, ]
  sort_cmd = " ".join(sort_cmd_list)
  DebugMsg(self, sort_cmd)
  exit_status = RunCommandFromString(sort_cmd, dpt=self.options.dpt)
  if (0 != exit_status): #{
    raise MyError("could not sort file with "
      "\"%s\": %i" % (" ".join(sort_cmd), exit_status))
  #} end if
#} end def
