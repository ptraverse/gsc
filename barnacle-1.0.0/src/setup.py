#! /usr/bin/env python
"""
setup.py

Created by Lucas Swanson
Copyright (c) 2012 Canada's Michael Smith Genome Sciences Centre. All rights reserved.
"""

# import standard modules
from optparse import OptionParser, OptionGroup
import os, re, sys, time, traceback

# ensure that the sys.path is set up to allow importing Barnacle modules
barnacle_dir = sys.path[0]
if ("" == barnacle_dir): #{
  barnacle_dir = os.getcwd()
#} end if
while (not os.path.isfile(os.path.join(barnacle_dir, "barnacle.pl"))): #{
  barnacle_dir = os.path.dirname(barnacle_dir)
  if ("" == barnacle_dir or "/" == barnacle_dir): #{
    print "cannot find BARNACLE directory"
    raise Exception
  #} end if
#} end while
barnacle_dir = os.path.expanduser(barnacle_dir)
barnacle_dir = os.path.abspath(barnacle_dir)
if (barnacle_dir not in sys.path): #{
  sys.path.insert(0, barnacle_dir)
#} end if

# import custom modules
from version import VERSION
from utils.log import GetLogPath, CloseLogFile
from utils.error import MyError
from utils.input import GetYesOrNoInput, GetStringInput
from utils.general import (SetupMainClass, TimeSpent, WriteCommand,
  CheckConfigCommands, GetCommand, SetupConfig)
from utils.messages import ErrMsg, LogMsg, DebugMsg, ExtremeDebugMsg
from utils.files_paths import (CheckFilePath, CheckDirPath, EnsureAbsPath,
  GetOutDir, FileBoxCls, CleanLine)
from utils.subprocesses import (RunCommandFromString, RunCommandFromList,
  STRING_OUT)

# CONSTANTS
ES_SUCCESS = 0
ES_RUN_ERR = 1
ES_OPT_ERR = 2
ES_PATH_ERR = 3
ES_EXCEPTION = 4
MSG_SUCCESS = "SUCCESS"
MSG_FAIL = "FAIL"

class BarnacleSetupCls: #{
  def __init__(self, options): #{
    options.setup_config = False
    SetupMainClass(self, options)
  #} end def

  def __del__(self): #{
    CloseLogFile(self)
  #} end def

  def Run(self): #{
    LogMsg(self, "Setting up Barnacle tool-suite...")
    start = time.time()
    self.SetupConfigFile()
    self.ReloadConfig()
    if (self.options.setup_annots): #{
      self.SetupAnnotations()
    else:
      LogMsg(self, "Skipping annotations setup.")
    #} end if
    self.SetupGapRealigner()
    LogMsg(self, "\nTime spent setting up Barnacle tool-suite: %s" %
      TimeSpent(start))
  #} end def

  def SetupConfigFile(self): #{
    #LogMsg(self, "In SetupConfigFile!")
    cfg_path = os.path.join(self.options.barnacle_src_dir, "barnacle.cfg")
    if (os.path.exists(cfg_path)): #{
      #LogMsg(self, "Config file barnacle.cfg already exists, would you "
      #    "like to overwrite it (y/n)? ", newline=False)
      # get user input (convert to upper-case first character)
      response = GetYesOrNoInput(self, "Config file barnacle.cfg already "
        "exists, would you like to overwrite it (y/n)? ")
      #LogMsg(self, "%s%s" % (prompt, response), write_to_screen=False)
      if ("N" == response): #{
        LogMsg(self, "Using existing barnacle.cfg file.")
        return
      #} end if
    else:
      LogMsg(self, "No barnacle.cfg file")
    #} end if
    start = time.time()
    LogMsg(self, "Setting up configuration file...")
    LogMsg(self, "-"*80)
    template_path = os.path.join(self.options.barnacle_src_dir,
        "template_barnacle.cfg")
    template_file = FileBoxCls(template_path, "r", skip_blank_lines=False)
    config_path = os.path.join(self.options.barnacle_src_dir,
        "barnacle.cfg")
    config_file = FileBoxCls(config_path, "w")
    section = None
    for template_line in template_file: #{
      if ("" == template_line or template_line.startswith("#")): #{
        LogMsg(self, "  %s" % template_line.replace("#", "Note: "))
        config_file.WriteLine(template_line)
        continue
      #} end if
      if (template_line.startswith("[")): #{
        config_file.WriteLine(template_line)
        section = template_line
        LogMsg(self, "  Setting up %s section..." % section)
        continue
      #} end if
      if ("[commands]" == section): #{
        config_line = self.GetCommandPath(template_line)
      elif ("[cluster]" == section):
        config_line = self.GetClusterValue(template_line)
      else:
        raise BarnacleSetupError("Unrecognized section: %s" % section)
      #} end if
      config_file.WriteLine(config_line)
    #} end for
    LogMsg(self,"-"*80)
    LogMsg(self, "Time spent setting up configuration file: %s" %
      TimeSpent(start))
  #} end def

  def GetCommandPath(self, template_line): #{
    cmd_name = template_line.split("=")[0]
    which_cmd = "which %s" % cmd_name
    if (cmd_name.endswith("_cluster")): #{
      which_cmd = "ssh %s \"%s\"" % (self.options.cluster_head,
          which_cmd.replace("_cluster",""))
    #} end if
    response = "N"
    cmd_path = CleanLine(RunCommandFromString(which_cmd,
      stdout=STRING_OUT)[0])
    if ("no %s in" % cmd_name.replace("_cluster","") in cmd_path): #{
      LogMsg(self, "    Could not find %s binary with \"which\":\n      %s" %
        (cmd_name, cmd_path))
    else:
      if (cmd_path.startswith("alias")): #{
        cmd_path = cmd_path.split()[-1]
      #} end if
      #LogMsg(self, "    Found %s binary at %s. Is this the correct path "
      #  "(y/n)? " % (cmd_name, cmd_path), newline=False)
      # get user input (convert to upper-case first character)
      response = GetYesOrNoInput(self, "    Found %s binary at %s. Is this "
        "the correct path (y/n)? " % (cmd_name, cmd_path))
      #LogMsg(self, response, write_to_screen=False)
    #} end if
    while ("N" == response): #{
      #LogMsg(self, "    Please enter the correct path for the %s binary: " %
      #  cmd_name, newline=False)
      cmd_path = GetStringInput(self, "    Please enter the correct path "
        "for the %s binary: " % cmd_name).split(";")[0]
      #LogMsg(self, cmd_path, write_to_screen=False)
      if (not os.path.isfile(cmd_path)): #{
        #LogMsg(self, "      There is no file at %s. Is this the correct path "
        #  "(y/n)? " % (cmd_path), newline=False)
        # get user input (convert to upper-case first character)
        response = GetYesOrNoInput(self, "      There is no file at %s. Is "
          "this the correct path (y/n)? " % (cmd_path))
        #LogMsg(self, response, write_to_screen=False)
      else:
        response = "Y"
      #} end if
    #} end if
    return "%s=%s" % (cmd_name, cmd_path)
  #} end def

  def GetClusterValue(self, template_line): #{
    field_name = template_line.split("=")[0]
    #LogMsg(self, "    Please enter a value for the %s cluster field: " %
    #  field_name, newline=False)
    field_value = GetStringInput(self, "    Please enter a value for the "
      "%s cluster field: " % field_name)
    #LogMsg(self, field_value, write_to_screen=False)
    return "%s=%s" % (field_name, field_value)
  #} end def

  def ReloadConfig(self): #{
    self.options.cfg = None
    SetupConfig(self)
  #} end def

  def SetupGapRealigner(self): #{
    start = time.time()
    CheckConfigCommands(self, "python")
    # get the appropriate paths
    gr_info = dict()
    self.PrepareToCompileGapRealigner(gr_info)
    LogMsg(self, "\nChecking gap_realigner...\n  Source: %s" %
        gr_info['src_path'])
    # if the binary does not exist, or the source has been modified more
    #  recently than the binary
    if (not os.path.exists(gr_info["bin_path"]) or
        NeedsUpdate(gr_info["bin_path"], gr_info["src_path"])): #{
      LogMsg(self, "  Compiling local binary...\n  %s" % gr_info['bin_path'])
      local_compile_cmd = gr_info['cmd_template'].replace("OUT",
        gr_info['bin_path'])
      # compile the code
      self.CompileGapRealigner(local_compile_cmd)
      LogMsg(self, "-"*80)
    else:
      LogMsg(self, "  Local binary up to date: %s" % gr_info['bin_path'])
    #} end if
    # if the binary does not exist, or the source has been modified more
    #  recently than the binary
    if (not os.path.exists(gr_info["cluster_bin_path"]) or
        NeedsUpdate(gr_info["cluster_bin_path"], gr_info["src_path"])): #{
      LogMsg(self, "  Compiling cluster binary...\n  %s" %
          gr_info['cluster_bin_path'])
      cluster_compile_cmd = "ssh %s \"%s\"" % (self.options.cluster_head,
        gr_info['cmd_template'].replace("OUT", gr_info['cluster_bin_path']))
      # compile the code for the cluster
      self.CompileGapRealigner(cluster_compile_cmd)
      LogMsg(self, "-"*80)
    else:
      LogMsg(self, "  Cluster binary up to date: %s" %
          gr_info['cluster_bin_path'])
    #} end if
    LogMsg(self, "Time spent compiling gap_realigner: %s" %
      TimeSpent(start))
  #} end def

  def PrepareToCompileGapRealigner(self, gr_info): #{
    gr_info['bin_path'] = os.path.join(
      self.options.barnacle_src_dir, "alignment_processing",
      "gap_realigner")
    gr_info['cluster_bin_path'] = ("%s_cluster" %
      gr_info['bin_path'])
    gr_info['src_path'] = ("%s.cpp" %
      gr_info['bin_path'])
    CheckFilePath(gr_info['src_path'], "gap_realigner source")
    gr_info['cmd_template'] = ("g++ -Wall -Werror -O3 -o OUT %s" %
      (gr_info['src_path']))
  #} end def

  def CompileGapRealigner(self, compile_command): #{
    LogMsg(self, compile_command)
    compile_status = RunCommandFromString(compile_command)
    if (0 > compile_status): #{
      raise BarnacleSetupError("Compile command was terminated by signal "
        "%i" % compile_status)
    elif (0 < compile_status):
      raise BarnacleSetupError("Error running compile command: %i" %
        compile_status)
    #} end if
  #} end def

  def SetupAnnotations(self): #{
    start = time.time()
    LogMsg(self, "\nSetting up annotation files...")
    LogMsg(self, "-"*80)
    parent_dir = os.path.dirname(self.options.barnacle_src_dir)
    annots_dir = os.path.join(parent_dir, "annotations")
    setup_annots_path = os.path.join(annots_dir, "setup_annotations.sh")
    # run the setup_annotations.sh script to download files from UCSC, etc.
    RunCommandFromString(setup_annots_path)
    # get the path to the create_gene_feature_coords.py script
    python = GetCommand(self, "python")
    create = os.path.join(self.options.barnacle_src_dir, "annotation",
        "create_gene_feature_coords.py")
    # get the list of gene feature coordinate files to create
    list_path = os.path.join(annots_dir, "gene_features_to_create.txt")
    for genes_file_name in FileBoxCls(list_path, "r"): #{
      in_path = os.path.join(annots_dir, genes_file_name)
      out_path = in_path.replace("txt", "exons.introns.std_chr.bed")
      if (NeedsUpdate(out_path, in_path)): #{
        command = [python, create, in_path, annots_dir]
        LogMsg(self, "")
        status = RunCommandFromList(command)
        CheckStatus(status, "create_gene_feature_coords.py")
      else:
        LogMsg(self, "\nGene feature coordinates file up to date:\n  %s" %
            out_path)
      #} end if
    #} end for
    LogMsg(self,"-"*80)
    LogMsg(self, "Time spent setting up annotation files: %s" %
      TimeSpent(start))
  #} end def
#} end class

def NeedsUpdate(dest_path, source_path): #{
  if (not os.path.exists(dest_path)): #{
    return True
  #} end if
  try: #{
    dest_time = os.path.getmtime(dest_path)
    source_time = os.path.getmtime(source_path)
  except OSError, e:
    raise BarnacleSetupError("Error checking file modification time: %s" % e)
  #} end try
  return (dest_time < source_time)
#} end def

def CheckStatus(status, descriptor): #{
  if (0 > status): #{
    raise BarnacleSetupError("%s command was terminated by signal %i" %
      (descriptor, -status))
  elif (0 < status):
    raise BarnacleSetupError("Error running %s command: %i" %
      (descriptor, status))
  #} end if
#} end def

#### EXCEPTION CLASSES ####
class BarnacleSetupError(MyError): #{
  pass
#} end class

def SetupOptionsParser(): #{
  description_string = ("Setup the Barnacle tool-suite: compile Gap "
    "Realigner binaries")
  args = [ "CLUSTER_HEAD", ]
  usage_string = "%prog " + " ".join(args) + " [ OPTIONS ]"
  parser = OptionParser(description=description_string,
                        version="%prog " + VERSION,
                        usage=usage_string)
  parser.num_args = len(args)
  parser.add_option("--no-annots",
                    action="store_false", dest="setup_annots",
                    help="Skipping setting up annotation files.")
  parser.add_option("-f", "--force",
                    action="store_true",
                    help="Force filtering to take place, even if the output "
                         "directory already exists.")
  parser.add_option("-d", "--debug",
                    action="store_true",
                    help="Print debug information while the program runs.")
  parser.set_defaults(setup_annots=True,
                      force=False,
                      debug=False)
  return parser
#} end def

def CheckPaths(options): #{
  #opts_good = True
  #if (not opts_good): #{
  #  ErrMsg("bad option") #TODO
  #  opts_good = False
  ##} end if
  #path_errors = list()
  #CheckFilePath(options.input_path, "INPUT", path_errors) #TODO
  ## get and check the output path
  #options.output_dir = GetOutDir(input_dir, "TASK_DESC") #TODO
  #if (opts_good and 0 == len(path_errors)): #{
  #  CheckDirPath(options.output_dir, "output", path_errors,
  #    create=True, replace=options.force)
  #  # get the log file name
  #  options.log_file_name = GetLogPath(options.input_path,
  #    "TASK_DESC", options.output_dir) #TODO
  ##} end if
  #if (0 < len(path_errors)): #{
  #  ErrMsg("Errors in input arguments:\n  %s" % "\n  ".join(path_errors))
  ##} end if
  ## the paths are good if there are no path errors and no conflicting options
  #return (opts_good and 0 == len(path_errors))
  options.log_file_name = os.path.join(barnacle_dir, "setup.log")
  return True
#} end def

def Main(): #{
  # get options and arguments
  parser = SetupOptionsParser()
  (options, args) = parser.parse_args()
  # if the right number of args was used
  if (parser.num_args == len(args)): #{
    options.cluster_head = args[0]
    if (CheckPaths(options)): #{
      try: #{
        main_class_object = BarnacleSetupCls(options)
        WriteCommand(main_class_object, sys.argv)
        main_class_object.Run()
      except (MyError), e:
        ErrMsg("ERROR while setting up Barnacle tool-suite:\n  %s" % e)
        return ES_RUN_ERR
      #} end try
    else:
      return ES_PATH_ERR
    #} end if
  else:
    parser.error("you must specify a cluster head to use to compile the "
      "gap_realigner tool for cluster use (CLUSTER_HEAD).")
    return ES_OPT_ERR
  #} end if
  return ES_SUCCESS
#} end def

if __name__ == '__main__': #{
  try: #{
    exit_status = Main()
  except Exception, e:
    traceback.print_exc()
    exit_status = ES_EXCEPTION
  #} end try
  if (ES_SUCCESS == exit_status): #{
    print MSG_SUCCESS
  else:
    print MSG_FAIL
  #} end if
  sys.exit(exit_status)
#} end if
