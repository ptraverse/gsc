#! /usr/bin/env python
"""
general.py

Created by Lucas Swanson
Copyright (c) 2012 Canada's Michael Smith Genome Sciences Centre. All rights reserved.
"""

# import standard modules
import itertools, math, os, re, string, sys, time, math
try: #{
  import ConfigParser
except ImportError:
  import configparser as ConfigParser
#} end try

# import custom modules
from log import OpenLogFile
from error import MyError
from messages import LogMsg, DebugMsg, ExtremeDebugMsg
from files_paths import EnsureAbsPath, CheckFilePath, PathError, FileBoxCls
from subprocesses import RunCommandFromString, STREAM_OUT

# constants
CFG_COMMANDS = "commands"
CFG_CLUSTER = "cluster"
CHR_ID_PATT = r"\A(chr)?([1-9][0-9]*|[XY]|MT?)\Z"
COORD_PATT = "%s:[0-9]+" % CHR_ID_PATT[:-2]
COORD_RANGE_PATT = "%s-[0-9]+" % COORD_PATT
LEFT  = "left"
RIGHT = "right"
SIDES = [LEFT, RIGHT]

def SetupMainClass(self, options, params=None,
    log_file=None, log_info=None): #{
  # set general functions
  self.options = options
  #self.debug = options.debug
  if (hasattr(options, "extreme_debug") and options.extreme_debug): #{
    options.debug = True
  #} end if
  opts_string = "{%s}" % ", ".join([
    r"'%s': %s" % (field, getattr(options, field)) for
    field in sorted(vars(options).keys())])
  if (None != params): #{
    opts_string = str(params)
    self.params = params
  #} end if
  SetupLogInfo(self, options, opts_string, log_file, log_info)
  # ensure that the BARNACLE_PATH environmental variable is set
  SetupEnviron(self, options)
  # assume that the config file should be read
  if (not hasattr(options, "setup_config") or options.setup_config): #{
    SetupConfig(self)
  #} end if
  ExtremeDebugMsg(self, "Finished SetupMainClass")
#} end def

def SetupLogInfo(self, options, opts_string, log_file, log_info): #{
  self.close_log_file = False
  write_options = True
  if (None != log_info): #{
    if ('file' in log_info): #{
      log_file = log_info['file']
    #} end if
    if ('write_opts' in log_info): #{
      write_options = log_info['write_opts']
    #} end if
  #} end if
  # if no log file is given
  if (None == log_file): #{
    # attempt to open log file
    if (hasattr(options, 'log_file_name')): #{
      self.log_file = OpenLogFile(options.log_file_name)
      self.close_log_file = True
    else:
      self.log_file = None
    #} end if
    # write options being used to log file
    if (write_options): #{
      if (options.debug): #{
        LogMsg(self, opts_string)
      elif (None != self.log_file): #{
        self.log_file.WriteLine(opts_string)
      #} end if
    #} end if
  else:
    # use the given log file
    self.log_file = log_file
  #} end if
  self.log_info = {'debug': options.debug, 'file': self.log_file,
    'write_opts': False, 'extreme_debug': False}
  if (hasattr(options, 'extreme_debug')): #{
    self.log_info['extreme_debug'] = options.extreme_debug
    if (self.log_info['extreme_debug']): #{
      self.log_info['debug'] = True
    #} end if
  #} end if
#} end def

def SetupEnviron(self, options): #{
  if ("BARNACLE_PATH" in os.environ): #{
    options.barnacle_src_dir = os.environ["BARNACLE_PATH"]
    return
  #} end if
  barnacle_dir = sys.path[0]
  DebugMsg(self, "Initially: %s" % barnacle_dir)
  if ("" == barnacle_dir): #{
    barnacle_dir = os.getcwd()
    DebugMsg(self, "Using current working dir: %s" % barnacle_dir)
  #} end if
  while (not os.path.isfile(os.path.join(barnacle_dir, "barnacle.pl"))): #{
    barnacle_dir = os.path.dirname(barnacle_dir)
    DebugMsg(self, "Up one dir: %s" % barnacle_dir)
    if ("" == barnacle_dir or "/" == barnacle_dir): #{
      raise MyError("cannot find BARNACLE directory")
    #} end if
  #} end while
  os.environ["BARNACLE_PATH"] = EnsureAbsPath(barnacle_dir)
  options.barnacle_src_dir = os.environ["BARNACLE_PATH"]
  DebugMsg(self, "BARNACLE_PATH set to %s" % barnacle_dir)
#} end def

def SetupConfig(self): #{
  if (hasattr(self.options, "cfg") and None != self.options.cfg): #{
    return
  #} end if
  DebugMsg(self, "Loading BARNACLE config file...")
  config_path = os.path.join(os.environ["BARNACLE_PATH"], "barnacle.cfg")
  CheckFilePath(config_path, "BARNACLE configuration")
  cfg = ConfigParser.SafeConfigParser()
  cfg.readfp(open(config_path))
  if (not cfg.has_section(CFG_COMMANDS)): #{
    raise MyError("config file must contain %s section\nPath: %s" %
      (CFG_COMMANDS, config_path))
  #} end if
  cfg.path = config_path
  self.options.cfg = cfg
#} end def

def ConfigHasCommand(self, command): #{
  return self.options.cfg.has_option(CFG_COMMANDS, command)
#} end def

# required is command, or list of commands, to check
def CheckConfigCommands(self, required, local=True, check_path=True,
cluster_head=None): #{
  if (isinstance(required, str)): #{
    required = [required]
  #} end if
  missing_commands = list()
  for command in required: #{
    if (not local): #{
      if (None != cluster_head): #{
        cluster_command = command + "_%s" % cluster_head
        try: #{
          CheckConfigCommands(self, cluster_command,
            local=True, check_path=check_path)
          continue
        except (CommandMissingError, PathError), e:
          DebugMsg(self, "Checking general cluster command for %s" % command)
        #} end try
      #} end if
      cluster_command = command + "_cluster"
      try: #{
        CheckConfigCommands(self, cluster_command,
          local=True, check_path=check_path)
        continue
      except (CommandMissingError, PathError), e:
        DebugMsg(self, "Using local command for %s" % command)
      #} end try
    #} end if
    if (ConfigHasCommand(self, command)): #{
      if (check_path): #{
        path = GetCommand(self, command, local=True)
        CheckFilePath(path, "%s command" % command)
      #} end if
    else:
      missing_commands.append(command)
    #} end if
  #} end for
  if (0 < len(missing_commands)): #{
    raise CommandMissingError("The following required commands are missing "
      "from the config file: %s\nConfig path: %s" %
      (", ".join(missing_commands), self.options.cfg.path))
  #} end if
#} end def

def GetCommand(self, command, local=True, cluster_head=None): #{
  if (not local): #{
    if (None != cluster_head): #{
      cluster_command = command + "_%s" % cluster_head
      if (ConfigHasCommand(self, cluster_command)): #{
        return self.options.cfg.get(CFG_COMMANDS, cluster_command)
      #} end if
    #} end if
    cluster_command = command + "_cluster"
    if (ConfigHasCommand(self, cluster_command)): #{
      return self.options.cfg.get(CFG_COMMANDS, cluster_command)
    #} end if
  #} end if
  return self.options.cfg.get(CFG_COMMANDS, command)
#} end def

def GetClusterValue(self, fieldname): #{
  if (self.options.cfg.has_option(CFG_CLUSTER, fieldname)): #{
    return self.options.cfg.get(CFG_CLUSTER, fieldname)
  #} end if
  raise ConfigError("configuration file has no \"%s\" value" % fieldname)
#} end def

def GetGroupID(member_id): #{
  group_pattern = r"^(\d+)[a-z]+"
  group_match   = re.search(group_pattern, member_id)
  if (None == group_match): #{
    return None
  #} end if
  return int(group_match.group(1))
#} end def

def IntOrNAString(data): #{
  if ("N/A" == data): #{
    return data
  else:
    return "%i" % data
  #} end if
#} end def

def ConvertSeconds(seconds): #{
  days    = int(math.floor(seconds / 60.0 / 60.0 / 24.0))
  hours   = int(math.floor(seconds / 60.0 / 60.0 % 24.0))
  minutes = int(math.floor(seconds / 60.0 % 60.0))
  seconds = seconds % 60.0
  time_parts = list()
  if (0 < days): #{
    time_parts.append("%id" % days)
  #} end if
  if (0 < hours): #{
    time_parts.append("%ih" % hours)
  #} end if
  if (0 < minutes): #{
    time_parts.append("%im" % minutes)
  #} end if
  time_parts.append("%is" % seconds)
  return " ".join(time_parts)
#} end def

def NormalizeGeneID(raw_gid): #{
  return raw_gid.replace(".non_coding","")
#} end def

def IsUTR(feature_name): #{
  return (feature_name.lower().endswith(".utr") or
    "utr" == feature_name.lower())
#} end def

def IsNotUTR(feature_name): #{
  return not IsUTR(feature_name)
#} end def

def IsNonCoding(feature_name): #{
  return feature_name.lower().endswith(".non_coding")
#} end def

def IsNotNonCoding(feature_name): #{
  return not IsNonCoding(feature_name)
#} end def

def IsNotUTRorNonCoding(feature_name): #{
  return (IsNotUTR(feature_name) and IsNotNonCoding(feature_name))
#} end def

def RemoveExonTypeLabel(feature_name): #{
  return feature_name.replace(".UTR","").replace(".non_coding","")
#} end def

def RemoveIndexingFlag(feature_name): #{
  return re.sub(r"X\d+X$", "", feature_name)
#} end def

#def CleanExonName(exon_name): #{
#  # remove the count
#  exon_name = re.sub(r"X\d+X$", "", exon_name)
#  # do NOT remove the part after the period, it is the exon type label
#  # (e.g. UTR or non_coding)
#  return exon_name.replace(" ", "_")
##} end def

def RunOverlapCode(target_coords_path, query_coords_path, output_path,
    get_size=False, get_closest=False, full=False, dpt=False, log_info=None): #{
  # create the overlap command
  overlap_code = os.path.join(os.environ["BARNACLE_PATH"],
    "utils", "overlapcoordinates_ultrafast.pl")
  overlap_cmd = "%s -coord %s -ref %s -us" % (overlap_code,
    target_coords_path, query_coords_path)
  if (get_size): #{
    overlap_cmd += " -size"
  #} end if
  if (get_closest): #{
    overlap_cmd += " -closest_dist"
  #} end if
  if (full): #{
    overlap_cmd += " -full"
  #} end if
  # setup the output file information
  fail_msg = "cannot open exon overlaps output file"
  if (STREAM_OUT == output_path): #{
    overlap_out_file = output_path
  else:
    overlap_out_file = {'path':output_path, 'mode':"w", 'fail_msg':fail_msg}
  #} end if
  DebugMsg(log_info, overlap_cmd)
  # call the overlap command
  return RunCommandFromString(overlap_cmd, stdout=overlap_out_file,
    dpt=dpt)
#} end def

def TimeSpent(start_time, calc_elapsed=True): #{
  if (calc_elapsed): #{
    elapsed = time.time() - start_time
  else:
    elapsed = start_time
  #} end if
  hours   = int(math.floor(elapsed / 60 / 60))
  minutes = int(math.floor(elapsed / 60 % 60))
  seconds = elapsed % 60
  time_parts = list()
  if (0 < hours): #{
    time_parts.append("%ih" % hours)
  #} end if
  if (0 < minutes): #{
    time_parts.append("%im" % minutes)
  #} end if
  time_parts.append("%.3fs" % seconds)
  return "".join(time_parts)
#} end def

# return the reverse complement of a DNA sequence
def ReverseComplement(s): #{
  t = string.maketrans("AaCcGgTt", "TtGgCcAa")
  return s.translate(t)[::-1]
#} end def

class OptsCls: #{
  def __init__(self): #{
    self.debug = True
  #} end def
#} end class

def GetOptions(module): #{
  (options, arguments) = module.SetupOptionsParser().parse_args([""])
  return options
#} end def

def CheckOptions(options, required_member_list): #{
  missing_members = list()
  for member in required_member_list: #{
    if (not hasattr(options, member)): #{
      missing_members.append(member)
    #} end if
  #} end for
  if (0 < len(missing_members)): #{
    raise MyError("The following required members are missing:\n  %s" %
      "\n  ".join(missing_members))
  #} end if
#} end def

def IsEmpty(list_to_check): #{
  if (0 == len(list_to_check)): #{
    return True
  #} end if
  if (1 == len(list_to_check)): #{
    if ("n/a"  == list_to_check[0].lower() or
        "none" == list_to_check[0].lower()):
      return True
    #} end if
  #} end if
  return False
#} end def

class ConfigError(MyError): #{
  pass
#} end class

class CommandMissingError(MyError): #{
  pass
#} end class

# check whether chromosome IDs in the given column of
# the given file use "chr"
def ShouldChromUseChr(column, path, file_description, log_info=None): #{
  fail_msg = "could not read %s coordinates file" % file_description
  file = FileBoxCls(path, "r", fail_msg)
  line = file.Head()
  chrom_value = line.split(" ")[column-1]
  if ("chr" in chrom_value): #{
    use_chr = True
  else:
    use_chr = False
  #} end if
  DebugMsg(log_info, "Line: \"%s\"\n" % line +
    "Chromosome Value: \"%s\" Use \"chr\": %s" % (chrom_value, use_chr))
  return use_chr
#} end def

def AddChr(target): #{
  if (target.startswith("chr")): #{
    return target
  elif (target.lower().startswith("chr")):
    return "chr%s" % target[3:]
  #} end if
  return "chr%s" % target
#} end if

def NormalizeChrID(raw_chr, use_chr=False): #{
  ouput_chr = raw_chr.upper().replace("CHR", "").replace("MT", "M")
  if (use_chr): #{
    ouput_chr = AddChr(ouput_chr)
  #} end if
  return ouput_chr
#} end def

def ChrToInt(chr): #{
  # make sure that the chromosome is just a single number or letter
  test_chr = NormalizeChrID(chr)
  letter_chrs = ["X", "Y", "M"]
  if (test_chr in letter_chrs): #{
    return -letter_chrs.index(test_chr)
  #} end if
  return int(test_chr)
#} end def

def NonStandardChr(chr): #{
  # make sure that the chromosome is just a single number or letter
  test_chr = NormalizeChrID(chr)
  return (None == re.search(CHR_ID_PATT, test_chr))
#} end def

def IsInteger(value): #{
  try: #{
    int_value = int(value)
  except ValueError:
    return False
  #} end try
  return True
#} end def

def WriteCommand(self, command_list): #{
  LogMsg(self, "COMMAND: %s" % " ".join(command_list), write_to_screen=False)
#} end def

def MaxMatchLength(query, target): #{
  max_match_len = 0
  # initialize last_row to a row of zeroes
  last_row = [0]*(len(target)+1)
  # initialize current_row to a single zero
  current_row = [0]
  #print "    %s" % " ".join(target)
  prev_char = " "
  for qend in range(len(query)): #{
    #print "%s %s" % (prev_char, " ".join(map(str,last_row)))
    prev_char = query[qend]
    for tend in range(len(target)): #{
      if (query[qend] == target[tend]): #{
        current_row.append(last_row[tend]+1)
        if (current_row[-1] > max_match_len): #{
          max_match_len = current_row[-1]
        #} end if
      else:
        current_row.append(0)
      #} end if
    #} end for
    last_row = current_row
    current_row = [0]
  #} end for
  #print "%s %s" % (query[qend], " ".join(map(str,last_row)))
  return max_match_len
#} end def

def IsHomopolymerSequence(sequence): #{
  first = sequence[0]
  for character in sequence[1:]: #{
    if (character != first): #{
      return False
    #} end if
  #} end for
  return True
#} end def

def PowerSet(iterable):
  "powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
  s = list(iterable)
  return itertools.chain.from_iterable(itertools.combinations(s, r) for
    r in range(len(s)+1))
#} end def

def OneMismatchRegexPattern(input_string): #{
  pattern_list = (input_string[:i]+"[ACGNT]".replace(c,'')+input_string[i+1:]
    for (i,c) in enumerate(input_string.upper()))
  return "|".join(pattern_list)
#} end def

def OtherSide(side): #{
  if (LEFT == side): #{
    return RIGHT
  elif (RIGHT == side):
    return LEFT
  else:
    raise MyError("unrecognized side: \"%s\"" % side)
  #} end if
#} end def

def StrJoin(delim, to_join): #{
  return delim.join(map(str, to_join))
#} end def

def IntFloor(value): #{
  return int(math.floor(value))
#} end def

def IntCeiling(value): #{
  return int(math.ceil(value))
#} end def

def Flatten(an_iterable): #{
  return itertools.chain.from_iterable(an_iterable)
#} end def

# from http://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#On-line_algorithm
def online_variance(data): #{
  n = 0
  mean = 0
  M2 = 0
  for x in data: #{
    n = n + 1
    delta = x - mean
    mean = mean + float(delta)/float(n)
    M2 = M2 + delta*(x - mean)
  #} end for
  variance_n = float(M2)/float(n)
  variance = float(M2)/float(n-1)
  return (variance, variance_n)
#} end def

# from http://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#On-line_algorithm
def online_std_dev(data): #{
  return math.sqrt(online_variance(data)[1])
#} end def
