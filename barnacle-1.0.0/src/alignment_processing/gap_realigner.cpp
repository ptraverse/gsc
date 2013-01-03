#include<iostream>
#include<fstream>
#include<sstream>
#include<string>
#include<cmath>
#include <iomanip>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
using namespace std;

string ver = "v1.11.0";

// Gap Realigner v1.11
// Tool for detecting chimeric events by re-aligning gaps in BLAT alignment to the contig and intronic sequence between blocks
//
// Created by Deniz Yorukoglu
// Edited by Lucas Swanson
// Copyright (c) 2012 Canada's Michael Smith Genome Sciences Centre. All rights reserved.

// 1.11
// Lucas: changed "TYPE" output field name to "TOPOLOGY"
//        increased precision of PID output

// 1.10
// Fixed a bug in the TandemCheck function

// 1.9
// Extends alignment region for duplications that are adjacent to the alignment neighbourhood boundary

// 1.8
// Added extraMetaData option
// Fixed a bug regarding reverse strand contigs that have a seperate block between the event and dup_event sequences

// 1.7
// Fixed alignment extension to avoid over-extension and incorrect non-tandem classification
// Fixed error in partial alignment threshold
// Fixed an inconsistency error between gap and duplicated regions.

// 1.6
// Fixed inversion block error

// 1.5
// Fixed output PID calculations
// Built seperate inversion search size threshold of 6
// Minimum Gap Threshold is changed to be as low as 3 (enforcing partial alignment minimum length as well)

// 1.4
// Detects and restores tandem spanning gaps
// Reports formations (Tandem Extension, Swapped ordering of gap and alignment, Tandem spanning gaps, Genomic Extension)

// 1.3
// Local alignment for duplications instead of BLAT

// 1.2
// Modified output according to output_spec_2

// 1.1
// Now it can detect duplications of sequences shorter than 20 bp
// And duplications found by partial alignment of gaps

// 1.0
// Tandem Duplication/Inv_Duplication -> Align the gap within the remaining contig
// Intronic_inversion -> align the gap to the the intron between adjacent blocks (reverse)

//These temporary file names will be added to output directory path and given mkstemp names.
string tempGapSeqFile = "_t_Gap";
string tempInvPSLOutputFile = "_t_PSL";
string tempHostNameFile = "_t_Host";
//temporary file for extracting sequences from genome
string tempExtractedStringFile = "_t_genExt";

//These files should not be deleted after program is killed, therefore only outputDir should be added and not mkstemp name
string junkBLAT_CMDout = "__temp_BLAT_cmd_out.txt";
string nonBasicEventLogFile = "nonBasicEvents.log";

string junkTwoBitExtractionErrorLogFileName = "__two_bit_conversion_error_log.txt";

//Command line parameters
//any gap >=threshold
int gapThreshold = 7;
double queryGapSimThreshold = 0.95;
double queryGapMinPartialMatchThresh = 0.4;
string extraMetaData = "";
bool debug(false);

//Config file parameters
string BlatSource = "";
string BlatExecPath = "";
string BlatHostName = "";
string BlatServer2BitDir = "";
string BlatClientPath = "";
string globalBlatOptions = "";
string globalBlatClientOptions = "";
string TwoBitGenPath = "";
string TwoBitToFaConverterPath = "";

string allGapRemovedQuery;
int relaxedVersionFLAG = 0;
int fullContigAlignFLAG = 0;
int tempFilesSetFLAG = 0;
int optionalGapPenalty = 1;
int CheckTSGFlag = 1;
int postPreFormCheckFLAG = 0;

//global data used by GapSearch and sub-procedures
int matches, mismatches, repmatches, junk, querySize, genStart, genEnd, genSize, blockCount, contigStart, contigEnd, qNumInsert;
string junks,queryName, chrID, contigDirection;
int *blockSizes, *qIndices, *tIndices;

int numNonBasicMatch = 0;

// Hardcoded Constant Parameters
const int negligibleShortTandemExtendingRegionThreshold = 2;
const double maxAllowedNonTandemDistanceRatio = 1.0;
const int negligibleGenomicDistanceThreshold = 5;
const int minGapLength_InversionSearch = 6;
const int absoluteMinGapLength = 3;

void DebugMsg(string msg) {
  if (debug) {
    cout << msg << "\n";
  }
}

void DebugMsg(stringstream msg) {
  DebugMsg(msg.str());
}

//UP(0)orDOWN(1) tells if the the sequence is upstream or downstream of the block w.r.t contig direction
int ConvertQueryPos2GenPos(int queryPos, int& trimLen, int UPorDOWN)
{
  DebugMsg("<ConvertQueryPos2GenPos>");
  trimLen = 0;
  if(contigDirection == "-") {
    //Coordinate is reversed for block direction
    queryPos = querySize - queryPos - 1;
    //Block end orientation is reversed
    UPorDOWN = 1 - UPorDOWN;
  }

  //Return the corresponding genome position of a query position
  for(int i=0; i<blockCount; i++) {
    if(queryPos < qIndices[i] + blockSizes[i]) {
      if(queryPos < qIndices[i]) {
        if(i==0 || UPorDOWN == 0) {
          // trimming up end
          trimLen = qIndices[i] - queryPos;
          DebugMsg("</ConvertQueryPos2GenPos>");
          return tIndices[i];
        } else {
          // trimming down end
          trimLen = queryPos - (qIndices[i-1] + blockSizes[i-1] - 1);
          DebugMsg("</ConvertQueryPos2GenPos>");
          return tIndices[i-1] + blockSizes[i-1] - 1;
        }
      }
      int offset = queryPos - qIndices[i];
      DebugMsg("</ConvertQueryPos2GenPos>");
      return tIndices[i] + offset;
    }
  }

  if(queryPos >= qIndices[blockCount-1] + blockSizes[blockCount-1] && queryPos < querySize) {
    //Always trim down end
    trimLen = queryPos - (qIndices[blockCount-1] + blockSizes[blockCount-1] - 1); 
    DebugMsg("</ConvertQueryPos2GenPos>");
    return tIndices[blockCount-1] + blockSizes[blockCount-1] - 1;
  }

  // if queryPos >= querySize
  cout << "ERROR Q: query position is larger than query size" << endl;
  exit(141);
}

int ConvertQueryPos2GenPos_Strict(int queryPos) //Strict version requires exact correspondence between query and genome indices, othewise return -1 (no trimming)
{
  if(contigDirection == "-") {
    //Coordinate is reversed for block direction
    queryPos = querySize - queryPos - 1;
  }

  //Return the corresponding genome position of a query position
  for(int i=0; i<blockCount; i++) {
    if(queryPos < qIndices[i] + blockSizes[i]) {
      if(queryPos < qIndices[i]) {
        return -1;
      }
      int offset = queryPos - qIndices[i];
      return tIndices[i] + offset;
    }
  }

  if(queryPos >= qIndices[blockCount-1] + blockSizes[blockCount-1] && queryPos < querySize) {
    return -1;
  }

  // if queryPos >= querySize
  cout << "ERROR Q: query position is larger than query size" << endl;
  exit(251);
}

int GetGenPos_LastBaseOfPreviousBlock(int queryPos)
{
  if(contigDirection == "-") {
    //Coordinate is reversed for block direction
    queryPos = querySize - queryPos - 1;
  }

  //Return genomic position of the last base of the previous block of query position
  //No block before
  if(queryPos < qIndices[0] + blockSizes[0]) {
    return -1;
  }

  for(int i=1; i<blockCount; i++) {
    if(queryPos < qIndices[i] + blockSizes[i]) {
      return tIndices[i-1] + blockSizes[i-1] - 1;
    }
  }

  return tIndices[blockCount-1] + blockSizes[blockCount-1];
}

int GetGenPos_FirstBaseOfFollowingBlock(int queryPos)
{
  if(contigDirection == "-") {
    //Coordinate is reversed for block direction
    queryPos = querySize - queryPos - 1;
  }

  //Return genomic position of the first base of following block of query position
  for(int i=0; i<blockCount; i++) {
    if(queryPos < qIndices[i]) {
      return tIndices[i];
    }
  }

  //No block after
  return -1;
}

int CheckIfCoordsInTheSameBlock(int coord1, int coord2)
{
  if(contigDirection == "-") {
    coord1 = querySize - coord1 - 1;
    coord2 = querySize - coord2 - 1;
  }

  if(coord2 < coord1) {
    int temp = coord1;
    coord1 = coord2;
    coord2 = temp;
  }

  for(int i=0; i<blockCount; i++) {
    if(coord1 < qIndices[i]) {
      return 0;
    } else {
      int blockEnd = qIndices[i] + blockSizes[i] -1;
      if(coord1 <= blockEnd) {
        if(coord2 <= blockEnd) {
          return 1;
        } else {
          return 0;
        }
      }
    }
  }
  return 0;
}

double calcPID(int qStart, int qEnd, int tStart, int tEnd, int qNumInsert, int qMatch, int qMisMatch, int qRepMatch)
{
  DebugMsg("<calcPID>");
  int milliBad = 0;

  int sizeMul = 1;
  int qAliSize, tAliSize, aliSize;
  int sizeDif;
  int insertFactor;
  int total;

  qAliSize = sizeMul * (qEnd - qStart);
  tAliSize = tEnd - tStart;
  aliSize = min(qAliSize, tAliSize);
  if (aliSize <= 0) {
    DebugMsg("</calcPID>");
    return 0;
  }
  sizeDif = qAliSize - tAliSize;
  if (sizeDif < 0) {
    sizeDif = 0;
  }

  insertFactor = qNumInsert;
  total = (sizeMul * (qMatch + qRepMatch + qMisMatch));
  if(total != 0) {
    int roundedValue = (int) round(3*log(1+sizeDif));
    milliBad = (1000 * (qMisMatch*sizeMul + insertFactor + roundedValue)) / total;
  }
  DebugMsg("</calcPID>");
  return 100.0 - (double)milliBad * 0.1 ;
}

double calcAF(int match, int mismatch, int repmatch, int size)
{
  return double(match + mismatch + repmatch) / double(size);
}

void CreateAllGapRemovedQuery(const string& query, int bCount, int* bSizes, int* bStarts, int qSize, string contigDirection)
{
  allGapRemovedQuery = query;

  if(contigDirection == "+") {
    for(int k=0; k<bStarts[0]; k++) {
      allGapRemovedQuery.at(k) = 'N';
    }

    for(int i=0; i<blockCount-1; i++) {
      int posGapStart = bStarts[i] + bSizes[i];
      if(posGapStart < bStarts[i+1]) {
        int posGapEnd = bStarts[i+1];
        for(int k=posGapStart; k<posGapEnd; k++) {
          allGapRemovedQuery.at(k) = 'N';
        }
      }
    }
    for(int k=bStarts[blockCount-1] + bSizes[blockCount-1]; k<qSize; k++) {
      allGapRemovedQuery.at(k) = 'N';
    }
  } else {
    //for first gap
    for(int k=qSize - bStarts[0]; k<qSize; k++) {
      allGapRemovedQuery.at(k) = 'N';
    }

    for(int i=0; i<blockCount-1; i++) {
      int posGapStart = bStarts[i] + bSizes[i];
      if(posGapStart < bStarts[i+1]) {
        int posGapEnd = bStarts[i+1];
        int removeStart = qSize - posGapEnd, removeEnd = qSize - posGapStart;
        for(int k=removeStart; k<removeEnd; k++) {
          allGapRemovedQuery.at(k) = 'N';
        }
      }
    }

    //for last gap
    for(int k=qSize - (bStarts[blockCount-1] + bSizes[blockCount-1]) - 1; k>=0; k--) {
      allGapRemovedQuery.at(k) = 'N';
    }
  }
}

void ParseBlockInfo(string bStr, int* bInfo, int bCount)
{
  int afterCommaIndex = 0;
  for(int i=0; i<bCount; i++) {
    int nextCommaPos = bStr.find(",", afterCommaIndex);
    string tempStr = bStr.substr(afterCommaIndex, nextCommaPos - afterCommaIndex);
    bInfo[i] = atoi(tempStr.c_str());
    afterCommaIndex = nextCommaPos + 1;
  }
}

bool doBasesMatch(char c1, char c2)
{
  //(c1 = c1 - 'a' + 'A')
  if(c1 >= 'a') {
    c1 -= 32;
  }
  if(c2 >= 'a') {
    c2 -= 32;
  }

  if(c1 == 'N' || c2 == 'N') {
    return 0;
  }

  if(c1==c2) {
    return 1;
  } else {
    return 0;
  }
}

bool doBasesRevCompMatch(char c1, char c2)
{
  if(c1 >= 'a') {
    c1 -= 32;
  }
  if(c2 >= 'a') {
    c2 -= 32;
  }

  if(c1 == 'N' || c2 == 'N') {
    return 0;
  }

  int sum = (int) c1 + (int) c2;

  //'A'+'T' or 'C'+'G'
  if(sum == 149 || sum == 138) {
    return 1;
  } else {
    return 0;
  }
}

int** alignMat;
int alignMatCapacity = 0;

int GetSequenceDifference(const string& seq1, const string& seq2, int orient)
{
  // Initialize or Resize alignment matrix if necessary
  int maxSeqLength = max(seq1.length(),seq2.length());
  if(alignMatCapacity < maxSeqLength + 1) {
    if(alignMatCapacity!=0) {
      for(int i=0; i<=alignMatCapacity; i++) {
        free(alignMat[i]);
      }
      free(alignMat);
    }

    alignMatCapacity = maxSeqLength * 2;
    alignMat = (int**) malloc ((alignMatCapacity+1)*sizeof(int*));
    for(int i=0; i<=alignMatCapacity; i++) {
      alignMat[i] = (int*) malloc ((alignMatCapacity+1)*sizeof(int));
    }
  }

  int size1 = seq1.length(), size2 = seq2.length();
  alignMat[0][0] = 0;
  for(int i=1; i<=maxSeqLength; i++) {
    alignMat[i][0] = alignMat[i-1][0] + 1;
    alignMat[0][i] = alignMat[0][i-1] + 1;
  }

  //Standard Needleman-Wunch - can add pruning to make it faster
  if(orient == 1) {
    for(int i=1; i<=size1; i++) {
      for(int j=1; j<=size2; j++) {
        alignMat[i][j] = min(min(alignMat[i-1][j-1] + !(doBasesMatch(seq1.at(i-1),seq2.at(j-1))), alignMat[i-1][j] + 1), alignMat[i][j-1] + 1);
      }
    }
    int positiveDirectionError = alignMat[size1][size2];
    return positiveDirectionError;
  } else {
    //orient==-1
    for(int i=1; i<=size1; i++) {
      for(int j=1; j<=size2; j++) {
        alignMat[i][j] = min(min(alignMat[i-1][j-1] + !(doBasesRevCompMatch(seq1.at(i-1),seq2.at(seq2.length()-1-(j-1)))), alignMat[i-1][j] + 1), alignMat[i][j-1] + 1);
      }
    }

    int negativeDirectionError = alignMat[size1][size2];
    return negativeDirectionError;
  }
}

bool CheckIfBlockExistsWithinInterval(int interStart, int interEnd)
{
  if(contigDirection == "-") {
    //Coordinate is reversed for block direction
    int newinterStart = querySize - interEnd - 1;
    int newinterEnd = querySize - interStart - 1;
    interStart = newinterStart;
    interEnd = newinterEnd;
  }

  for(int i=0; i<blockCount; i++) {
    if(qIndices[i] >= interStart && qIndices[i] + blockSizes[i] - 1 <= interEnd) {
      return 1;
    }
  }
  return 0;
}

// Assigns ideal gap indices as well but it might have minor shifting problems until double sided tandem check is introduced
bool isDuplicationTandemByLocalAlign(char gap2contigDirection, int gap_Size, int dupLengthForGap, int dupLengthForGAlign, int gap_StartInd, int gGapStartPos, int gAlignStartPos, const string& query, int& ideal_gapStartInd_offset, int& ideal_gapEndInd_offset, string& ctg_extFormationFlag)
{
  DebugMsg("<isDuplicationTandemByLocalAlign>");
  //overblock extension check (if there is a seperate block between the two intervals, do not apply tandem extension)
  ctg_extFormationFlag = "";

  int gAlignEndPos = gAlignStartPos + dupLengthForGAlign - 1;

  int gOrgStartPos = gap_StartInd + gGapStartPos;
  int gOrgEndPos = gOrgStartPos + dupLengthForGap - 1;

  int insIntervalStart, insIntervalEnd;
  //gap is on the right
  if(gAlignStartPos < gOrgStartPos) {
    insIntervalStart = gAlignEndPos + 1;
    insIntervalEnd = gOrgStartPos - 1;
  //gap is on the left
  } else {
    insIntervalStart = gOrgEndPos + 1;
    insIntervalEnd = gAlignStartPos - 1;
  }

  if(CheckIfBlockExistsWithinInterval(insIntervalStart, insIntervalEnd)) {
    DebugMsg("  Block in interval\n</isDuplicationTandemByLocalAlign>");
    return 0;
  }

  ideal_gapStartInd_offset = 0;
  ideal_gapEndInd_offset = 0;

  if(gap2contigDirection == '+') {
    int distanceBetween;
    //gap is on the right
    if(gAlignStartPos < gOrgStartPos) {
      distanceBetween = gOrgStartPos - gAlignEndPos - 1;
    //gap is on the left
    } else {
      distanceBetween = gAlignStartPos - gOrgEndPos - 1;
    }

    //Too long gap to be considered as tandem - this threshold might be modified
    if(distanceBetween > dupLengthForGap) {
      DebugMsg("  Too Long\n</isDuplicationTandemByLocalAlign>");
      return 0;
    }

    // Enough short gap to be considered tandem, don't need to check further
    if(distanceBetween <= gap_Size * (1 - queryGapSimThreshold)) {
      //extension would cause index trouble, no reason to verify tandem extending region and modify ideal gap positions
      if(gAlignEndPos+1 < 0 || gAlignEndPos+distanceBetween >= (int) query.length() || gAlignStartPos-distanceBetween < 0 || gAlignStartPos - 1 >= (int) query.length()) {
        ctg_extFormationFlag = "TAN(Sh&Ir)";
        DebugMsg("  TAN(Sh&Ir)\n</isDuplicationTandemByLocalAlign>");
        return 1;
      }

      // Extension sequence too short, so ideal positions assigned as extended without looking at the sequences.
      if(distanceBetween <= negligibleShortTandemExtendingRegionThreshold) {
        //gap is on the right - ideal gap should be extended to gap.
        if(gAlignStartPos < gOrgStartPos) {
          ideal_gapStartInd_offset = -distanceBetween;
        } else {
          ideal_gapEndInd_offset = distanceBetween;
        }

        ctg_extFormationFlag = "TAN(Sh&Neg)";
        DebugMsg("  TAN(Sh&Neg)\n</isDuplicationTandemByLocalAlign>");
        return 1;
      }

      string Region1 = query.substr(gAlignEndPos+1,distanceBetween);
      string Region2 = query.substr(gAlignStartPos-distanceBetween,distanceBetween);

      int sequenceDifference = GetSequenceDifference(Region1, Region2, 1);

      //Relaxed difference threshold within extending region, but very strict within overall gap
      if(sequenceDifference <= (int) (3 * distanceBetween * (1.0 - queryGapSimThreshold))) {
        //Extending sequences can be assumed to be similar - Modifying ideal positions
        //gap is on the right - ideal gap should be extended to gap.
        if(gAlignStartPos < gOrgStartPos) {
          ideal_gapStartInd_offset = -distanceBetween;
        } else {
          ideal_gapEndInd_offset = distanceBetween;
        }
        ctg_extFormationFlag = "TAN(Sh&Mod)";
      } else {
        ctg_extFormationFlag = "TAN(Sh&NMod)";
      }

      DebugMsg("  Tandem\n</isDuplicationTandemByLocalAlign>");
      //Since distance is short compared to overall gap, reports as tandem in all cases. (But ideal positions aren't modified in all cases)
      return 1;
    }

    //For every case in between, check if the regions to make them tandem are similar
    //Coordinates of these regions do not depend on whether the gap is on the right or left, they just switch places.
    if(gAlignEndPos+1 < 0 || gAlignEndPos+distanceBetween >= (int) query.length() || gAlignStartPos-distanceBetween < 0 || gAlignStartPos - 1 >= (int) query.length()) {
      DebugMsg("  Not Tandem\n</isDuplicationTandemByLocalAlign>");
      return 0;
    }

    string Region1 = query.substr(gAlignEndPos+1,distanceBetween);
    string Region2 = query.substr(gAlignStartPos-distanceBetween,distanceBetween);
    int sequenceDifference = GetSequenceDifference(Region1, Region2, 1);
    if(sequenceDifference <= (int) (3 * distanceBetween * (1.0 - queryGapSimThreshold))) {
      //gap is on the right - ideal gap should be extended to gap.
      if(gAlignStartPos < gOrgStartPos) {
        ideal_gapStartInd_offset = -distanceBetween;
      } else {
        ideal_gapEndInd_offset = distanceBetween;
      }
      ctg_extFormationFlag = "TAN(Mod)";
    } else {
      ctg_extFormationFlag = "TAN(NMod)";
    }

    if(sequenceDifference <= max(negligibleShortTandemExtendingRegionThreshold, (int) (gap_Size * (1.0 - queryGapSimThreshold)))) {
      DebugMsg("  Tandem\n</isDuplicationTandemByLocalAlign>");
      return 1;
    } else {
      //resetting offsets to 0, since they won't be extended
      ideal_gapStartInd_offset = 0;
      ideal_gapEndInd_offset = 0;
      DebugMsg("  Not Tandem\n</isDuplicationTandemByLocalAlign>");
      return 0;
    }
  } else if(gap2contigDirection == '-') {
    int distanceBetween;
    //gap is on the right
    if(gAlignStartPos < gOrgStartPos) {
      distanceBetween = gOrgStartPos - gAlignEndPos - 1;
    //gap is on the left
    } else {
      distanceBetween = gAlignStartPos - gOrgEndPos - 1;
    }

    //This is twices as large compared to forward, since in reverse tandem connecting regions are adjacent
    if(distanceBetween > 2 * dupLengthForGap) {
      DebugMsg("  Not Tandem\n</isDuplicationTandemByLocalAlign>");
      return 0;
    }
    // Again twice as large from palindromic property
    if(distanceBetween <= 2 * gap_Size * (1 - queryGapSimThreshold)) {
      if(distanceBetween <= 2 * negligibleShortTandemExtendingRegionThreshold) {
        if(gAlignEndPos < gOrgStartPos) {
          ideal_gapStartInd_offset = -(distanceBetween/2);
        } else {
          ideal_gapEndInd_offset = distanceBetween/2;
        }

        ctg_extFormationFlag = "TAN(Sh&Neg)";
        DebugMsg("  TAN(Sh&Neg)\n</isDuplicationTandemByLocalAlign>");
        return 1;
      }

      string Region1, Region2;
      //gap is on the right
      if(gAlignEndPos < gOrgStartPos) {
        Region1 = query.substr(gAlignEndPos+1, distanceBetween/2);
        Region2 = query.substr(gOrgStartPos - distanceBetween/2, distanceBetween/2);
      } else {
        Region1 = query.substr(gOrgEndPos+1, distanceBetween/2);
        Region2 = query.substr(gAlignStartPos - distanceBetween/2, distanceBetween/2);
      }

      int sequenceDifference = GetSequenceDifference(Region1, Region2, -1);

      //Relaxed difference threshold within extending region, but very strict within overall gap
      if(sequenceDifference <= (int) (3 * distanceBetween * (1.0 - queryGapSimThreshold))) {
        if(gAlignEndPos < gOrgStartPos) {
          ideal_gapStartInd_offset = -(distanceBetween/2);
        } else {
          ideal_gapEndInd_offset = distanceBetween/2;
        }
        ctg_extFormationFlag = "TAN(Sh&Mod)";
      } else {
        ctg_extFormationFlag = "TAN(Sh&NMod)";
      }

      DebugMsg("  Tandem\n</isDuplicationTandemByLocalAlign>");
      return 1;
    }
    //For every case in between, check if the regions to make them tandem are similar
    //If distance between is odd, the letter in the exact middle isn't counted
    string Region1, Region2;
    //gap is on the right
    if(gAlignEndPos < gOrgStartPos) {
      Region1 = query.substr(gAlignEndPos+1, distanceBetween/2);
      Region2 = query.substr(gOrgStartPos - distanceBetween/2, distanceBetween/2);
    } else {
      Region1 = query.substr(gOrgEndPos+1, distanceBetween/2);
      Region2 = query.substr(gAlignStartPos - distanceBetween/2, distanceBetween/2);
    }

    int sequenceDifference = GetSequenceDifference(Region1, Region2, -1);

    //Relaxed difference threshold within extending region, but very strict within overall gap
    if(sequenceDifference <= (int) (3 * distanceBetween * (1.0 - queryGapSimThreshold))) {
      if(gAlignEndPos < gOrgStartPos) {
        ideal_gapStartInd_offset = -(distanceBetween/2);
      } else {
        ideal_gapEndInd_offset = distanceBetween/2;
      }
      ctg_extFormationFlag = "TAN(Mod)";
    } else {
      ctg_extFormationFlag = "TAN(NMod)";
    }

    if(sequenceDifference <= max(negligibleShortTandemExtendingRegionThreshold , (int) (gap_Size * (1.0 - queryGapSimThreshold)) )) {
      DebugMsg("  Tandem\n</isDuplicationTandemByLocalAlign>");
      return 1;
    } else {
      ///resetting offsets to 0, since they won't be extended
      ideal_gapStartInd_offset = 0;
      ideal_gapEndInd_offset = 0;
      DebugMsg("  Not Tandem (seq too different)\n</isDuplicationTandemByLocalAlign>");
      return 0;
    }
  } else {
    cout << "ERROR: unknown direction character - Probably caused by faulty PSL line" << endl;
    exit(138);
  }
}

int CheckTandemSpanningGap(const string& Query, int origGapStartInd, int origGapEndInd, int idealGapStartInd, int idealGapEndInd, int idealGAlignStartInd, int idealGAlignEndInd, int& idealGapStartInd_offset, int& idealGapEndInd_offset)
{
  DebugMsg("<CheckTandemSpanningGap>");
  idealGapStartInd_offset = 0;
  idealGapEndInd_offset = 0;

  // The sequences to be aligned are always on both side of gapDup and length is the distance between gapDup and gAlignDup
  int dupLengthForGap = idealGapEndInd - idealGapStartInd + 1;
  int distanceBetween;

  //gap is on the left side (compare the sequunce between gap & gAlign with before gap)
  if(idealGapEndInd < idealGAlignStartInd) {
    //In this case duplication within the gap should be adjacent to upstream border
    if(idealGapStartInd != origGapStartInd) {
      DebugMsg("</CheckTandemSpanningGap>");
      return 0;
    }

    distanceBetween = idealGAlignStartInd - idealGapEndInd - 1;
    //Too long gap to be considered as tandem - this threshold might be modified
    if(distanceBetween > dupLengthForGap) {
      DebugMsg("</CheckTandemSpanningGap>");
      return 0;
    }

    //there isn't enough space before to extend before gap
    if(idealGapStartInd < distanceBetween) {
      DebugMsg("</CheckTandemSpanningGap>");
      return 0;
    }

    string Region1 = Query.substr(idealGapStartInd - distanceBetween, distanceBetween);
    string Region2 = Query.substr(idealGapEndInd + 1, distanceBetween);

    int sequenceDifference = GetSequenceDifference(Region1, Region2, 1);

    if(sequenceDifference <= (int) (distanceBetween * (1.0 - queryGapSimThreshold)) ) {
      idealGapStartInd_offset = -distanceBetween;
      DebugMsg("</CheckTandemSpanningGap>");
      return 1;
    } else {
      DebugMsg("</CheckTandemSpanningGap>");
      return 0;
    }
  // gap is on the right side (compare the sequence
  } else if(idealGAlignEndInd < idealGapStartInd) {
    //In this case duplication within the gap should be adjacent to downstream border
    if(idealGapEndInd != origGapEndInd) {
      DebugMsg("</CheckTandemSpanningGap>");
      return 0;
    }

    distanceBetween = idealGapStartInd - idealGAlignEndInd - 1;

    //Too long gap to be considered as tandem - this threshold might be modified
    if(distanceBetween > dupLengthForGap) {
      DebugMsg("</CheckTandemSpanningGap>");
      return 0;
    }

    if(idealGapEndInd + distanceBetween >= (int) Query.length()) {
      DebugMsg("</CheckTandemSpanningGap>");
      return 0;
    }

    string Region1 = Query.substr(idealGapStartInd - distanceBetween, distanceBetween);
    string Region2 = Query.substr(idealGapEndInd + 1, distanceBetween);

    int sequenceDifference = GetSequenceDifference(Region1, Region2, 1);

    if(sequenceDifference <= (int) (distanceBetween * (1.0 - queryGapSimThreshold)) ) {
      idealGapEndInd_offset = distanceBetween;
      DebugMsg("</CheckTandemSpanningGap>");
      return 1;
    } else {
      DebugMsg("</CheckTandemSpanningGap>");
      return 0;
    }
  } else {
    cout << "ERROR: gapDup and gAlignDup sequences shouldn't be overlapping" << endl;
    exit(140);
  }
  DebugMsg("  reached end\n</CheckTandemSpanningGap>");
}

void PrintGapDuplicationSplitAlignmentGenomeCoordinates(int g2genomeIndexStart,int g2genomeIndexEnd, ofstream& fout, int bCount, const int* bSizes, const int* tInd)
{
  DebugMsg("<PrintGapDuplicationSplitAlignmentGenomeCoordinates>");
  // If indices are given reverse, switch them
  if(g2genomeIndexStart > g2genomeIndexEnd) {
    int temp = g2genomeIndexStart;
    g2genomeIndexStart = g2genomeIndexEnd;
    g2genomeIndexEnd = temp;
  }

  for(int i=0; i<bCount; i++) {
    if(g2genomeIndexStart < tInd[i] && tInd[i] + bSizes[i] - 1 < g2genomeIndexEnd) {
      fout << tInd[i]+1 << "-" << tInd[i] + bSizes[i] << ",";
    } else {
      if(g2genomeIndexStart >= tInd[i] && tInd[i] + bSizes[i] - 1 >= g2genomeIndexEnd) {
        fout << g2genomeIndexStart + 1 << "-" << g2genomeIndexEnd + 1;
        DebugMsg("</PrintGapDuplicationSplitAlignmentGenomeCoordinates>");
        return;
      } else if(g2genomeIndexStart >= tInd[i] && g2genomeIndexStart <= tInd[i] + bSizes[i] - 1) {
        fout << g2genomeIndexStart + 1 << "-" << tInd[i] + bSizes[i] << ",";
      } else if(g2genomeIndexEnd >= tInd[i] && g2genomeIndexEnd <= tInd[i] + bSizes[i] - 1) {
        fout << tInd[i]+1 << "-" << g2genomeIndexEnd + 1;
        DebugMsg("</PrintGapDuplicationSplitAlignmentGenomeCoordinates>");
        return;
      }
    }
  }
}

int getBaseMatchScore(char c1, char c2)
{
  if(doBasesMatch(c1,c2)) {
    return 1;
  } else {
    return -1;
  }
}

int getBaseRevCompMatchScore(char c1, char c2)
{
  if(doBasesRevCompMatch(c1,c2)) {
    return 1;
  } else {
    return -1;
  }
}

int DetermineMaxAlignLenAndStartPos(int ctgLen, int gapLen, int gapStart, int& before_startInd, int& after_startInd, int& before_partLen, int& after_partLen)
{
  int maxAllowed = (int) ((double) gapLen * ( 1 + maxAllowedNonTandemDistanceRatio));
  before_startInd = gapStart - maxAllowed;
  before_partLen = maxAllowed;
  after_startInd = gapStart + gapLen;
  after_partLen = maxAllowed;

  int beforeGapLen = gapStart;
  if(beforeGapLen < maxAllowed) {
    before_startInd = 0;
    before_partLen = beforeGapLen;
  }

  int afterGapLen = ctgLen - (gapStart + gapLen);
  if(afterGapLen < maxAllowed) {
    after_partLen = afterGapLen;
  }

  int maxAftBef = max(beforeGapLen, afterGapLen);
  if(maxAftBef <= maxAllowed) {
    return maxAftBef;
  }
  return maxAllowed;
}

void GetMaxScoreAndInterval(int** localAlignMat, const string& gap_seq, const string& contig_seq, int ctgStart, char gap_dir, int endX, int endY, int& startX, int& startY)
{
  int posX = endX, posY=endY;
  while(localAlignMat[posX][posY]!=0)
  {
    int isSame;
    if(gap_dir == '+') {
      isSame = doBasesMatch(gap_seq.at(posX-1),contig_seq.at(ctgStart+posY-1));
    } else {
      isSame = doBasesRevCompMatch(gap_seq.at(gap_seq.length()-posX),contig_seq.at(ctgStart+posY-1));
    }

    //1 is diagonal, 2 is up, 3 is side
    int bestPreDir = 1;
    int bestPreScore = localAlignMat[posX-1][posY-1];

    //bactracking towards up and left can only happen if bases don't match
    if(isSame != 1) {
      int sideScore = localAlignMat[posX][posY-1];
      if(sideScore > bestPreScore) {
        bestPreScore = sideScore;
        bestPreDir = 3;
      }

      int upScore = localAlignMat[posX-1][posY];
      if(upScore > bestPreScore) {
        bestPreScore = upScore;
        bestPreDir = 2;
      }
    }

    startX = posX;
    startY = posY;

    switch(bestPreDir)
    {
      case 1:
        posX = posX-1;
        posY = posY-1;
        break;
      case 2:
        posX = posX-1;
        break;
      case 3:
        posY = posY-1;
        break;
      default:
        break;
    }
  }
}

int LocalAlignQueryWithinFullContig(const string& gapSeq, const string& contigSeq, int gap_StartInd, char& direction, int& withinGap_DupStartPos, int& withinGap_DupEndPos, int& gAlign_DupStartPos, int& gAlign_DupEndPos, int&  dup_MatchScore)
{
  DebugMsg("<LocalAlignQueryWithinFullContig>");
  int gapSeqLen = gapSeq.length(), contigSeqLen = contigSeq.length();
  int gapAlignmentScoreThreshold;// = gapThreshold;
  int gapPenalty = optionalGapPenalty;

  if(relaxedVersionFLAG == 1) {
    gapAlignmentScoreThreshold = (int) (gapSeqLen * queryGapMinPartialMatchThresh);
  } else {
    gapAlignmentScoreThreshold = (int) (gapSeqLen * queryGapSimThreshold);
  }

  if(gapAlignmentScoreThreshold < gapThreshold) {
    gapAlignmentScoreThreshold = gapThreshold;
  }

  int** locAlignMat = (int**) malloc ((gapSeqLen+1)*sizeof(int*));
  for(int i=0; i<=gapSeqLen; i++) {
    locAlignMat[i] = (int*) malloc ((contigSeqLen+1)*sizeof(int));
  }

  //These are common for all 4 alignments since they aren't edited at all
  for(int i=0; i<=gapSeqLen; i++) {
    locAlignMat[i][0] = 0;
  }
  for(int j=1; j<=contigSeqLen; j++) {
    locAlignMat[0][j] = 0;
  }

  int forwardScore = gapAlignmentScoreThreshold-1, endX=-1, endY=-1;
  for(int i=1; i<= gapSeqLen; i++) {
    for(int j=1; j<= contigSeqLen; j++) {
      locAlignMat[i][j] = max(0,max(max( locAlignMat[i-1][j-1] + getBaseMatchScore(gapSeq.at(i-1),contigSeq.at(j-1)), locAlignMat[i-1][j] - gapPenalty), locAlignMat[i][j-1] - gapPenalty));
      if(locAlignMat[i][j]>forwardScore) {
        forwardScore = locAlignMat[i][j];
        endX = i;
        endY = j;
      }
    }
  }

  int forwardGapIntStart=-1, forwardCtgIntStart=-1, forwardCtgIntEnd=-1, forwardGapIntEnd=-1;
  if(forwardScore >= gapAlignmentScoreThreshold) {
    int startX=-1, startY=-1;
    GetMaxScoreAndInterval(locAlignMat, gapSeq, contigSeq, 0, '+', endX, endY, startX, startY);

    forwardGapIntStart = startX - 1;
    forwardGapIntEnd = endX -1;

    forwardCtgIntStart = startY - 1;
    forwardCtgIntEnd = endY -1;
  }

  int reverseScore = gapAlignmentScoreThreshold-1; endX=-1; endY=-1;
  for(int i=1; i<=gapSeqLen; i++) {
    for(int j=1; j<=contigSeqLen; j++) {
      locAlignMat[i][j] = max(0,max(max( locAlignMat[i-1][j-1] + getBaseRevCompMatchScore(gapSeq.at(gapSeqLen-i),contigSeq.at(j-1)), locAlignMat[i-1][j] - gapPenalty), locAlignMat[i][j-1] - gapPenalty) );
      if(locAlignMat[i][j]>reverseScore) {
        reverseScore = locAlignMat[i][j];
        endX = i;
        endY = j;
      }
    }
  }

  int reverseGapIntStart=-1, reverseGapIntEnd=-1, reverseCtgIntStart=-1, reverseCtgIntEnd=-1;
  if(reverseScore >= gapAlignmentScoreThreshold) {
    int startX=-1, startY=-1;
    GetMaxScoreAndInterval(locAlignMat, gapSeq, contigSeq, 0, '-', endX, endY, startX, startY);

    reverseGapIntStart = gapSeqLen - endX;
    reverseGapIntEnd = gapSeqLen - startX;

    reverseCtgIntStart = startY - 1;
    reverseCtgIntEnd = endY -1;
  }

  for(int i=0; i<=gapSeqLen; i++) {
    free(locAlignMat[i]);
  }
  free(locAlignMat);

  // 1 for BF (Before forward), 2 for AF (after forward), 3 for BR (before reverse) , 4 for AR (after reverse)
  int maxScore = forwardScore; int maxEvent = 1;
  int maxScoringLen = forwardGapIntEnd - forwardGapIntStart + 1;

  if(reverseScore > maxScore) {
    maxScore = reverseScore;
    maxEvent = 2;
    maxScoringLen = reverseGapIntEnd - reverseGapIntStart + 1;
  }

  if(maxScore < gapAlignmentScoreThreshold) {
    DebugMsg("  Score Too Low\n</LocalAlignQueryWithinFullContig>");
    return 0;
  }

  ///Check Similarity Score w.r.t. interval length
  double simScorePart = (double) maxScore / maxScoringLen;

  if(simScorePart < queryGapSimThreshold) {
    DebugMsg("  Similarity Too Low\n</LocalAlignQueryWithinFullContig>");
    return 0;
  }

  dup_MatchScore = maxScore;
  if(maxEvent == 1) {
    withinGap_DupStartPos = forwardGapIntStart;
    withinGap_DupEndPos = forwardGapIntEnd;
    gAlign_DupStartPos = forwardCtgIntStart;
    gAlign_DupEndPos = forwardCtgIntEnd;
    direction = '+';
    DebugMsg("  Max Event == 1\n</LocalAlignQueryWithinFullContig>");
    return 1;
  } else {
    withinGap_DupStartPos = reverseGapIntStart;
    withinGap_DupEndPos = reverseGapIntEnd;
    gAlign_DupStartPos = reverseCtgIntStart;
    gAlign_DupEndPos = reverseCtgIntEnd;
    direction = '-';
    DebugMsg("  Max Event != 1\n</LocalAlignQueryWithinFullContig>");
    return 1;
  }
}

int RealignGapWithNewSeq(const string& gapSeq, const string& contigSeq, int newStartPos, int newEndPos, char direction, int prevMaxScore, int& ret_gapIntStart, int& ret_gapIntEnd, int& ret_ctgIntStart, int &ret_ctgIntEnd) //return matchScore
{
  int gapSeqLen = gapSeq.length(), contigPartLen = newEndPos - newStartPos + 1;
  int gapPenalty = optionalGapPenalty;
  int** locAlignMat = (int**) malloc ((gapSeqLen+1)*sizeof(int*));
  for(int i=0; i<=gapSeqLen; i++) {
    locAlignMat[i] = (int*) malloc ((contigPartLen+1)*sizeof(int));
  }
  for(int i=0; i<=gapSeqLen; i++) {
    locAlignMat[i][0] = 0;
  }
  for(int j=1; j<=contigPartLen; j++) {
    locAlignMat[0][j] = 0;
  }
  int bestScore = prevMaxScore, endX=-1, endY=-1;;
  if(direction == '+') {
    for(int i=1; i<=gapSeqLen; i++) {
      for(int j=1; j<=contigPartLen; j++) {
        locAlignMat[i][j] = max(0,max(max( locAlignMat[i-1][j-1] + getBaseMatchScore(gapSeq.at(i-1),contigSeq.at(newStartPos+j-1)), locAlignMat[i-1][j] - gapPenalty), locAlignMat[i][j-1] - gapPenalty) );
        if(locAlignMat[i][j] > bestScore) {
          bestScore = locAlignMat[i][j];
          endX = i;
          endY = j;
        }
      }
    }
    if(bestScore > prevMaxScore) {
      int startX=-1, startY=-1;
      GetMaxScoreAndInterval(locAlignMat, gapSeq, contigSeq, newStartPos, '+', endX, endY, startX, startY);

      ret_gapIntStart = startX - 1;
      ret_gapIntEnd = endX -1;

      ret_ctgIntStart = newStartPos + startY - 1;
      ret_ctgIntEnd = newStartPos + endY -1;
    }
  } else {
    for(int i=1; i<= gapSeqLen; i++) {
      for(int j=1; j<= contigPartLen; j++) {
        locAlignMat[i][j] = max(0,max(max( locAlignMat[i-1][j-1] + getBaseRevCompMatchScore(gapSeq.at(gapSeqLen-i),contigSeq.at(newStartPos+j-1)), locAlignMat[i-1][j] - gapPenalty), locAlignMat[i][j-1] - gapPenalty) );
        if(locAlignMat[i][j]>bestScore) {
            bestScore = locAlignMat[i][j];
            endX = i;
            endY = j;
        }
      }
    }
    if(bestScore > prevMaxScore) {
        int startX=-1, startY=-1;
        GetMaxScoreAndInterval(locAlignMat, gapSeq, contigSeq, newStartPos, '-', endX, endY, startX, startY);

        ret_gapIntStart = gapSeqLen - endX;
        ret_gapIntEnd = gapSeqLen - startX;

        ret_ctgIntStart = newStartPos + startY - 1;
        ret_ctgIntEnd = newStartPos + endY -1;
    }
  }

  for(int i=0; i<=gapSeqLen; i++) {
    free(locAlignMat[i]);
  }
  free(locAlignMat);

  return bestScore;
}

int LocalAlignQueryWithinNeighboringArea(const string& gapSeq, const string& contigSeq, int gap_StartInd, char& direction, int& withinGap_DupStartPos, int& withinGap_DupEndPos, int& gAlign_DupStartPos, int& gAlign_DupEndPos, int&  dup_MatchScore)
{
  DebugMsg("<LocalAlignQueryWithinNeighboringArea>");
  int gapSeqLen = gapSeq.length(), contigSeqLen = contigSeq.length();

  int gapAlignmentScoreThreshold;// = gapThreshold;
  int gapPenalty = optionalGapPenalty;

  if(relaxedVersionFLAG == 1) {
    gapAlignmentScoreThreshold = (int) (gapSeqLen * queryGapMinPartialMatchThresh);
  } else {
    gapAlignmentScoreThreshold = (int) (gapSeqLen * queryGapSimThreshold);
  }

  if(gapAlignmentScoreThreshold < gapThreshold) {
    gapAlignmentScoreThreshold = gapThreshold;
  }

  int beforeStart, afterStart, beforePartLen, afterPartLen;
  int contigPartLen = DetermineMaxAlignLenAndStartPos(contigSeqLen, gapSeqLen, gap_StartInd, beforeStart, afterStart, beforePartLen, afterPartLen);

  int** locAlignMat = (int**) malloc ((gapSeqLen+1)*sizeof(int*));
  for(int i=0; i<=gapSeqLen; i++) {
    locAlignMat[i] = (int*) malloc ((contigPartLen+1)*sizeof(int));
  }
  //These are common for all 4 alignments since they aren't edited at all
  for(int i=0; i<=gapSeqLen; i++) {
    locAlignMat[i][0] = 0;
  }
  for(int j=1; j<=contigPartLen; j++) {
    locAlignMat[0][j] = 0;
  }

  int BFscore = gapAlignmentScoreThreshold-1, endX=-1, endY=-1;

  for(int i=1; i<= gapSeqLen; i++) {
    for(int j=1; j<= beforePartLen; j++) {
      locAlignMat[i][j] = max(0,max(max( locAlignMat[i-1][j-1] + getBaseMatchScore(gapSeq.at(i-1),contigSeq.at(beforeStart+j-1)), locAlignMat[i-1][j] - gapPenalty), locAlignMat[i][j-1] - gapPenalty) );
      if(locAlignMat[i][j]>BFscore) {
        BFscore = locAlignMat[i][j];
        endX = i;
        endY = j;
      }
    }
  }

  int BFgapIntStart=-1, BFctgIntStart=-1, BFctgIntEnd=-1, BFgapIntEnd=-1;
  if(BFscore >= gapAlignmentScoreThreshold) {
    int startX=-1, startY=-1;
    GetMaxScoreAndInterval(locAlignMat, gapSeq, contigSeq, beforeStart, '+', endX, endY, startX, startY);

    BFgapIntStart = startX - 1;
    BFgapIntEnd = endX -1;

    BFctgIntStart = beforeStart + startY - 1;
    BFctgIntEnd = beforeStart + endY -1;
  }

  int AFscore = gapAlignmentScoreThreshold-1;
  for(int i=1; i<= gapSeqLen; i++) {
    for(int j=1; j<= afterPartLen; j++) {
      locAlignMat[i][j] = max(0,max(max( locAlignMat[i-1][j-1] + getBaseMatchScore(gapSeq.at(i-1),contigSeq.at(afterStart+j-1)), locAlignMat[i-1][j] - gapPenalty), locAlignMat[i][j-1] - gapPenalty) );
      if(locAlignMat[i][j]>AFscore) {
        AFscore = locAlignMat[i][j];
        endX = i;
        endY = j;
      }
    }
  }

  int AFgapIntStart = -1, AFgapIntEnd = -1, AFctgIntStart = -1, AFctgIntEnd = -1;
  if(AFscore >= gapAlignmentScoreThreshold) {
    int startX = -1, startY = -1;
    GetMaxScoreAndInterval(locAlignMat, gapSeq, contigSeq, afterStart, '+', endX, endY, startX, startY);
    AFgapIntStart = startX - 1;
    AFgapIntEnd = endX - 1;
    AFctgIntStart = afterStart + startY - 1;
    AFctgIntEnd = afterStart + endY -1;
  }

  int BRscore = gapAlignmentScoreThreshold-1;
  for(int i=1; i<= gapSeqLen; i++) {
    for(int j=1; j<= beforePartLen; j++) {
      locAlignMat[i][j] = max(0,max(max( locAlignMat[i-1][j-1] + getBaseRevCompMatchScore(gapSeq.at(gapSeqLen-i),contigSeq.at(beforeStart+j-1)), locAlignMat[i-1][j] - gapPenalty), locAlignMat[i][j-1] - gapPenalty) );
      if(locAlignMat[i][j]>BRscore) {
        BRscore = locAlignMat[i][j];
        endX = i;
        endY = j;
      }
    }
  }

  int BRgapIntStart=-1, BRgapIntEnd=-1, BRctgIntStart=-1, BRctgIntEnd=-1;
  if(BRscore >= gapAlignmentScoreThreshold) {
    int startX=-1, startY=-1;
    GetMaxScoreAndInterval(locAlignMat, gapSeq, contigSeq, beforeStart, '-', endX, endY, startX, startY);

    BRgapIntStart = gapSeqLen - endX;
    BRgapIntEnd = gapSeqLen - startX;

    BRctgIntStart = beforeStart + startY - 1;
    BRctgIntEnd = beforeStart + endY -1;
  }

  int ARscore = gapAlignmentScoreThreshold-1;
  for(int i=1; i<= gapSeqLen; i++) {
    for(int j=1; j<= afterPartLen; j++) {
      locAlignMat[i][j] = max(0,max(max( locAlignMat[i-1][j-1] + getBaseRevCompMatchScore(gapSeq.at(gapSeqLen-i),contigSeq.at(afterStart+j-1)), locAlignMat[i-1][j] - gapPenalty), locAlignMat[i][j-1] - gapPenalty) );
      if(locAlignMat[i][j]>ARscore) {
        ARscore = locAlignMat[i][j];
        endX = i;
        endY = j;
      }
    }
  }

  int ARgapIntStart=-1, ARgapIntEnd=-1, ARctgIntStart=-1, ARctgIntEnd=-1;
  if(ARscore >= gapAlignmentScoreThreshold) {
    int startX=-1, startY=-1;
    GetMaxScoreAndInterval(locAlignMat, gapSeq, contigSeq, afterStart, '-', endX, endY, startX, startY);

    ARgapIntStart = gapSeqLen - endX;
    ARgapIntEnd = gapSeqLen - startX;

    ARctgIntStart = afterStart + startY - 1;
    ARctgIntEnd = afterStart + endY -1;
  }

  // Later, instead of deallocating matrix for each gap, reallocate the matrix if a larger gap is going to be searched
  for(int i=0; i<=gapSeqLen; i++) {
    free(locAlignMat[i]);
  }
  free(locAlignMat);

  // 1 for BF (Before forward), 2 for AF (after forward), 3 for BR (before reverse) , 4 for AR (after reverse)
  int maxScore = BFscore; int maxEvent = 1;
  int maxScoringLen = BFgapIntEnd - BFgapIntStart + 1;

  if(AFscore > maxScore) {
    maxScore = AFscore;
    maxEvent = 2;
    maxScoringLen = AFgapIntEnd - AFgapIntStart + 1;
  }
  if(BRscore > maxScore) {
    maxScore = BRscore;
    maxEvent = 3;
    maxScoringLen = BRgapIntEnd - BRgapIntStart + 1;
  }
  if(ARscore > maxScore) {
    maxScore = ARscore;
    maxEvent = 4;
    maxScoringLen = ARgapIntEnd - ARgapIntStart + 1;
  }

  if(maxScore < gapAlignmentScoreThreshold) {
    DebugMsg("  Score Too Low\n</LocalAlignQueryWithinNeighboringArea>");
    return 0;
  }

  ///Check Similarity Score w.r.t. interval length
  double simScorePart = (double) maxScore / maxScoringLen;
  if(simScorePart < queryGapSimThreshold) {
    DebugMsg("  Similarity Too Low\n</LocalAlignQueryWithinNeighboringArea>");
    return 0;
  }

  dup_MatchScore = maxScore;
  switch(maxEvent)
  {
    case 1:
      withinGap_DupStartPos = BFgapIntStart;
      withinGap_DupEndPos = BFgapIntEnd;
      gAlign_DupStartPos = BFctgIntStart;
      gAlign_DupEndPos = BFctgIntEnd;
      direction = '+';
            break;
    case 2:
      withinGap_DupStartPos = AFgapIntStart;
      withinGap_DupEndPos = AFgapIntEnd;
      gAlign_DupStartPos = AFctgIntStart;
      gAlign_DupEndPos = AFctgIntEnd;
      direction = '+';
            break;
    case 3:
      withinGap_DupStartPos = BRgapIntStart;
      withinGap_DupEndPos = BRgapIntEnd;
      gAlign_DupStartPos = BRctgIntStart;
      gAlign_DupEndPos = BRctgIntEnd;
      direction = '-';
             break;
    case 4:
      withinGap_DupStartPos = ARgapIntStart;
      withinGap_DupEndPos = ARgapIntEnd;
      gAlign_DupStartPos = ARctgIntStart;
      gAlign_DupEndPos = ARctgIntEnd;
      direction = '-';
      break;
    default:
      DebugMsg("</LocalAlignQueryWithinNeighboringArea>");
      return 0;
  }

  if(maxEvent == 1 || maxEvent == 3) {
          ///Check if beginning of the alignment falls to the starting edge of the BFinterval
    if(gAlign_DupStartPos == beforeStart) {
              ///Define the new extended alignment region and re-align (noting direction of alignment)
              int newStartPos = max(0, beforeStart - gapSeqLen);
              int newEndPos = max(gAlign_DupEndPos, beforeStart + gapSeqLen - 1);

              int ret_gapIntStart, ret_gapIntEnd, ret_ctgIntStart, ret_ctgIntEnd;
              int score = RealignGapWithNewSeq(gapSeq, contigSeq, newStartPos, newEndPos, direction, dup_MatchScore, ret_gapIntStart, ret_gapIntEnd, ret_ctgIntStart, ret_ctgIntEnd);

              ///Store new alignment values here
              if(score > dup_MatchScore) {
                  dup_MatchScore = score;
                  withinGap_DupStartPos = ret_gapIntStart;
                  withinGap_DupEndPos = ret_gapIntEnd;
                  gAlign_DupStartPos = ret_ctgIntStart;
                  gAlign_DupEndPos = ret_ctgIntEnd;
              }
    }
  } else if(maxEvent == 2 || maxEvent == 4) {
          ///Check if the end of the alignment falls to the ending edge of the AFinterval
    int afterEnd = afterStart +  afterPartLen - 1;
    if(gAlign_DupEndPos == afterEnd) {
              ///Define the new extended alignment region and re-align
              int newEndPos = min(int(contigSeq.length()-1), afterEnd + gapSeqLen);
              int newStartPos = min(gAlign_DupStartPos, afterEnd - gapSeqLen + 1);

              int ret_gapIntStart, ret_gapIntEnd, ret_ctgIntStart, ret_ctgIntEnd;
              int score = RealignGapWithNewSeq(gapSeq, contigSeq, newStartPos, newEndPos, direction, dup_MatchScore, ret_gapIntStart, ret_gapIntEnd, ret_ctgIntStart, ret_ctgIntEnd);

              ///Store new alignment values here
              if(score > dup_MatchScore) {
                  dup_MatchScore = score;
                  withinGap_DupStartPos = ret_gapIntStart;
                  withinGap_DupEndPos = ret_gapIntEnd;
                  gAlign_DupStartPos = ret_ctgIntStart;
                  gAlign_DupEndPos = ret_ctgIntEnd;
              }
          }
  }

  DebugMsg("</LocalAlignQueryWithinNeighboringArea>");
  if(maxEvent >=1 && maxEvent <= 4) {
    return 1;
  } else {
      return 0;
    }
}

int OneSidedLocalAlignment(const string& genString, const string& ctgString, int& newEnd)
{
  if(genString.length() != ctgString. length()) {
    cout << "ERROR: Lengths of compared sequences for genomic extension isn't equal. (genomeString: " << genString.length() << "  contigString: " << ctgString.length() << ")" << endl;
    exit(136);
  }

  // Initialize or Resize alignment matrix if necessary
  int seqLength = genString.length();
  if(alignMatCapacity < seqLength + 1) {
    if(alignMatCapacity!=0) {
      for(int i=0; i<=alignMatCapacity; i++) {
        free(alignMat[i]);
      }
      free(alignMat);
    }

    alignMatCapacity = seqLength * 2;
    alignMat = (int**) malloc ((alignMatCapacity+1)*sizeof(int*));
    for(int i=0; i<=alignMatCapacity; i++) {
      alignMat[i] = (int*) malloc ((alignMatCapacity+1)*sizeof(int));
    }
  }

  alignMat[0][0] = 0;
  for(int i=1; i<=seqLength; i++) {
    alignMat[i][0] = alignMat[i-1][0] - 1;
    alignMat[0][i] = alignMat[0][i-1] - 1;
  }

  int l1GapPenalty = 1;
  //Modified Needleman-Wunch (Scores are like in Smith-waterman, however resetting scores at 0 isn't allowed, values can go negative)
  int maxScoreIndGen=0, maxScoreIndCtg=0;
  for(int i=1; i<=seqLength; i++) {
    for(int j=1; j<=seqLength; j++) {
      alignMat[i][j] = max(max(alignMat[i-1][j-1] + getBaseMatchScore(genString.at(i-1),ctgString.at(j-1)), alignMat[i-1][j] - l1GapPenalty), alignMat[i][j-1] - l1GapPenalty);
      if(alignMat[i][j] > alignMat[maxScoreIndGen][maxScoreIndCtg]) {
        maxScoreIndGen = i;
        maxScoreIndCtg = j;
      }
    }
  }

  //The index in the genome is accepted as the extension border
  newEnd = maxScoreIndGen;
  return alignMat[maxScoreIndGen][maxScoreIndCtg];
}

void Add_OutputDir_And_Mkstemp_To_TempFiles(const string& outputFileName)
{
  tempFilesSetFLAG = 1;

  size_t last_slash_pos = outputFileName.find_last_of("/");
  if(last_slash_pos != string::npos) {
    string outputDir(outputFileName.substr(0,last_slash_pos+1));
    tempGapSeqFile = outputDir + tempGapSeqFile;
    tempInvPSLOutputFile = outputDir + tempInvPSLOutputFile;
    tempHostNameFile = outputDir + tempHostNameFile;
    tempExtractedStringFile = outputDir + tempExtractedStringFile;
    junkBLAT_CMDout = outputDir + junkBLAT_CMDout;
    nonBasicEventLogFile = outputDir + nonBasicEventLogFile;
  }

  tempGapSeqFile += ".XXXXXX";
  tempInvPSLOutputFile += ".XXXXXX";
  tempHostNameFile += ".XXXXXX";
  tempExtractedStringFile += ".XXXXXX";

  char temp[200];
  strcpy(temp, tempGapSeqFile.c_str());
  mkstemp(temp);
  string tempgapstr(temp);
  tempGapSeqFile = tempgapstr;

  strcpy(temp, tempInvPSLOutputFile.c_str());
  mkstemp(temp);
  string tempInvPSLstr(temp);
  tempInvPSLOutputFile = tempInvPSLstr;

  strcpy(temp, tempExtractedStringFile.c_str());
  mkstemp(temp);
  string tempExt(temp);
  tempExtractedStringFile = tempExt;

}

void CleanTempFiles()
{
  string removeCall = "rm " + tempInvPSLOutputFile;
  system(removeCall.c_str());
  removeCall = "rm " + tempGapSeqFile;
  system(removeCall.c_str());
  removeCall = "rm " + tempExtractedStringFile;
  system(removeCall.c_str());
}

void AssignTemporaryHostNameFile()
{
  char temp[200];
  strcpy(temp, tempHostNameFile.c_str());
  mkstemp(temp);
  string tempHostNameStr(temp);
  tempHostNameFile = tempHostNameStr;
}

int CountNsInGap(const string& gap_seq)
{
  int Ncount = 0, gap_len = gap_seq.length();
  for(int i=0; i<gap_len; i++) {
    if(gap_seq.at(i) == 'N' || gap_seq.at(i) == 'n') {
      Ncount++;
    }
  }
  return Ncount;
}

bool DoesGapOnlyAppearInContig(int blockIndex, int distThreshold)
{
  if(blockIndex == 0 || blockIndex == blockCount) {
    return 0;
  }
  if(tIndices[blockIndex] - (tIndices[blockIndex-1] + blockSizes[blockIndex-1]) <= distThreshold) {
    return 1;
  }
  return 0;
}

int CheckIfIntervalIsUnaligned(int insStart, int insEnd)
{
  DebugMsg("<CheckIfIntervalIsUnaligned>");
  int negligibleSideOverlap = (int) ( double(insEnd - insStart +1) * (1 - queryGapSimThreshold) / 2.0) + 1;

  if(insEnd < qIndices[0] || insStart >= qIndices[blockCount-1] + blockSizes[blockCount-1]) {
    DebugMsg("</CheckIfIntervalIsUnaligned>");
    return 1;
  }

  for(int i=1; i<blockCount; i++) {
    if(qIndices[i-1] + blockSizes[i-1] - negligibleSideOverlap > insStart) {
      DebugMsg("</CheckIfIntervalIsUnaligned>");
      return 0;
    }

    if(insEnd < qIndices[i] && insStart >= qIndices[i-1] + blockSizes[i-1]) {
      DebugMsg("</CheckIfIntervalIsUnaligned>");
      return 1;
    } else if(insEnd - negligibleSideOverlap < qIndices[i] && insStart + negligibleSideOverlap >= qIndices[i-1] + blockSizes[i-1]) {
      DebugMsg("</CheckIfIntervalIsUnaligned>");
      return 2;
    }
  }

  DebugMsg("</CheckIfIntervalIsUnaligned>");
  return 0;
}

string ExtractFrom2Bit(int start, int end) //start and end positions are inclusive - corrdinates should be modified for twoBit2Fa calls.
{
  stringstream twoBitExtractionCall;
  twoBitExtractionCall << TwoBitToFaConverterPath << " " << TwoBitGenPath << ":" << chrID << ":" << start << "-" << end+1 << " " << tempExtractedStringFile;
  system((twoBitExtractionCall.str()).c_str());

  ifstream finExtStr(tempExtractedStringFile.c_str());
  string extStr;
  //discard the first line
  getline(finExtStr, extStr);
  getline(finExtStr, extStr);
  return extStr;
}

char compChar(char ch)
{
  switch(ch)
  {
    case 'a':
      return 't';
    case 't':
      return 'a';
    case 'g':
      return 'c';
    case 'c':
      return 'g';
    case 'A':
      return 'T';
    case 'T':
      return 'A';
    case 'G':
      return 'C';
    case 'C':
      return 'G';
    case 'N':
      return 'N';
    case 'n':
      return 'n';
    default:
    {
      printf("ERROR: wrong character in sequence input: %c\n",ch);
      exit(135);
    }
  }
}

string TakeRevCompOfSequence(const string& origString)
{
  string revString = origString;
  int len = revString.length();
  for(int i=0; i<len; i++) {
    revString.at(i) = compChar(origString.at(len - i -1));
  }
  return revString;
}

// 0-> N/A, 1-> only Post, 2->only Pre, 3->post+pre, 4->post+novel+pre 5->post+novel 6->pre+novel
int CheckPostPreFormation(int gapDupStart, int gapDupEnd, int gAlignStart, int gAlignEnd, int origGapStartInd, int origGapEndInd, const string& fullCtgSeq, int& ret_postLength, int& ret_novelLength, int& ret_preLength)
{
  DebugMsg("<CheckPostPreFormation>");
  int preFlag = 0, novelFlag = 0, postFlag = 0;
  if(gapDupEnd < gAlignStart) {
    //Post interval is from gapdup to the gap's right boundary
    int postStart = gapDupEnd+1, postEnd = origGapEndInd;
    //pre interval is from right after original gap to the gAlignDup
    int preStart = origGapEndInd+1, preEnd = gAlignStart - 1;

    int postInterval = postEnd - postStart + 1, preInterval = preEnd - preStart + 1;

    if(preInterval >= negligibleShortTandemExtendingRegionThreshold && CheckIfCoordsInTheSameBlock(preStart, gAlignStart)) {
      preFlag = 1;
      ret_preLength = preInterval;
    }
    //Strict version requires an exact correspondence of indices, otherwise it will return -1
    int dupGenomeEnd =  ConvertQueryPos2GenPos_Strict(gAlignEnd);

    if(dupGenomeEnd != -1 && postInterval >= negligibleShortTandemExtendingRegionThreshold) {
      int PostExtractStart = dupGenomeEnd + 1, PostExtractEnd = dupGenomeEnd + postInterval;

      string postString = ExtractFrom2Bit(PostExtractStart, PostExtractEnd);

      if(contigDirection == "-") {
        postString = TakeRevCompOfSequence(postString);
      }

      string ctgPostString = fullCtgSeq.substr(postStart,postInterval);

      //Here the alignment between the genomeString and contigString should be done in a 1-sided local alignment fashion.
      //First bases should start together, scoring is same as local, for backtracking the highest score can be taken from any cell in the alignment.
      //However, resetting score isn't allowed, thus alignment score can go negative. This is to prevent any prefix shifts, which would indicate that Post-Pre formation won't apply.
      int newPostEndingPos = 0;
      int postAlignmentScore = OneSidedLocalAlignment(postString, ctgPostString, newPostEndingPos);

      //To count as post region, 1- should cover at least postInterval * partialMatchThesh , 2- should have a very high alignment score w.r.t. aligned region
      if(postAlignmentScore >= newPostEndingPos * queryGapSimThreshold && postAlignmentScore >= postInterval * queryGapMinPartialMatchThresh) {
        postFlag = 1;
        ret_postLength = postInterval;
        int novelIntervalDistance = postInterval - newPostEndingPos;
        if(novelIntervalDistance >=negligibleShortTandemExtendingRegionThreshold) {
          novelFlag = 1;
          ret_novelLength = novelIntervalDistance;
          // update the length of postInterval after novelInterval is separated.
          ret_postLength -= ret_novelLength;
        }
      }
    }
  } else {
    //Pre interval is from gap's left boundary up to gapdup
    int preStart = origGapStartInd, preEnd = gapDupStart - 1;
    //Post interval is from after the gAlignDup, until the gap's beginning
    int postStart = gAlignEnd + 1, postEnd = origGapStartInd -1;
    int postInterval = postEnd - postStart + 1, preInterval = preEnd - preStart + 1;

    if(postInterval >= negligibleShortTandemExtendingRegionThreshold && CheckIfCoordsInTheSameBlock(postEnd, gAlignEnd)) {
      postFlag = 1;
      ret_postLength = postInterval;
    }
    //Strict version requires an exact correspondence of indices, otherwise it will return -1
    int dupGenomeStart =  ConvertQueryPos2GenPos_Strict(gAlignStart);
    if(dupGenomeStart != -1 && preInterval >= negligibleShortTandemExtendingRegionThreshold) {
      int PreExtractStart = dupGenomeStart - preInterval, PreExtractEnd = dupGenomeStart - 1;
      string preString = ExtractFrom2Bit(PreExtractStart, PreExtractEnd);

      //Here the genome string is reversed if the direction is forward, since the string from contig is reversed in all cases for 1-side local alignment (so it's like a double reverse for '-' strand).
      if(contigDirection == "+") {
        preString = TakeRevCompOfSequence(preString);
      }

      string ctgPreString = fullCtgSeq.substr(preStart,preInterval);
      //Reverse complement should be taken for the contig sequence, since 1-sided local alignment aligns left side together and lets right end loose.
      ctgPreString = TakeRevCompOfSequence(ctgPreString);

      //Whereas in this case we would like to have right side aligned together.
      int newReversePreEndingPos = 0;
      int preAlignmentScore = OneSidedLocalAlignment(preString, ctgPreString, newReversePreEndingPos);

      //To count as pre region, 1- should cover at least preInterval * partialMatchThesh , 2- should have a very high alignment score w.r.t. aligned region
      //Thus relaxed within overall extension, however, very strict within the pre region within extension
      if(preAlignmentScore >= newReversePreEndingPos * queryGapSimThreshold && preAlignmentScore >= postInterval * queryGapMinPartialMatchThresh) {
        preFlag = 1;
        ret_preLength = preInterval;
        int novelIntervalDistance = preInterval - newReversePreEndingPos;
        if(novelIntervalDistance >=negligibleShortTandemExtendingRegionThreshold) {
          novelFlag = 1;
          ret_novelLength = novelIntervalDistance;
          // update the length of preInterval after novelInterval is separated.
          ret_preLength -= ret_novelLength;
        }
      }
    }
  }

  if(preFlag + postFlag == 0) {
    DebugMsg("  No post or pre\n</CheckPostPreFormation>");
    return 0;
  }

  if(preFlag + postFlag == 2) {
    if(novelFlag == 0) {
      DebugMsg("  PostPre 3\n</CheckPostPreFormation>");
      return 3;
    } else {
      DebugMsg("  PostPre 4\n</CheckPostPreFormation>");
      return 4;
    }
  }

  if(preFlag == 1) {
    if(novelFlag == 0) {
      DebugMsg("  PostPre 2\n</CheckPostPreFormation>");
      return 2;
    } else {
      DebugMsg("  PostPre 6\n</CheckPostPreFormation>");
      return 6;
    }
  //postFlag = 1
  } else {
    if(novelFlag == 0) {
      DebugMsg("  PostPre 1\n</CheckPostPreFormation>");
      return 1;
    } else {
      DebugMsg("  PostPre 5\n</CheckPostPreFormation>");
      return 5;
    }
  }
}

bool CheckPolyAorT(const string& gapStr)
{
  char firstChar = gapStr.at(0);
  int len = gapStr.length();

  if(firstChar == 't' || firstChar == 'T') {
    for(int i=1; i<len; i++) {
      char curChar = gapStr.at(i);
      if(!(curChar == 't' || curChar == 'T')) {
        return 0;
      }
    }
    return 1;
  }

  if(firstChar == 'a' || firstChar == 'A') {
    for(int i=1; i<len; i++) {
      char curChar = gapStr.at(i);
      if(!(curChar == 'a' || curChar == 'A')) {
        return 0;
      }
    }
    return 1;
  }

  return 0;
}

bool CheckSubsetPolyAorT(const string& gapStr, int startPos, int endPos) //start and end positions are both inclusive
{
  char firstChar = gapStr.at(startPos);
  int len = endPos - startPos + 1;

  if(firstChar == 't' || firstChar == 'T') {
    for(int i=1; i<len; i++) {
      char curChar = gapStr.at(startPos+i);
      if(!(curChar == 't' || curChar == 'T')) {
        return 0;
      }
    }
    return 1;
  }

  if(firstChar == 'a' || firstChar == 'A') {
    for(int i=1; i<len; i++) {
      char curChar = gapStr.at(startPos+i);
      if(!(curChar == 'a' || curChar == 'A')) {
        return 0;
      }
    }
    return 1;
  }

  return 0;
}



int GapSearch(stringstream& str, const string& fullQuery, const string& outputFN, const string& genesList, const int BLAT_port_code)
{
  DebugMsg("<GapSearch>");
  // Number of gap events found. This is the return value of GapSearch.
  int numFoundMatch = 0;

  // Loading necesseary values from PSL line to global variables
  str >> matches >> mismatches >> repmatches >> junk >> qNumInsert >> junk >> junk >> junk >> contigDirection >> queryName >> querySize >> contigStart >> contigEnd >> chrID >> genSize >> genStart >> genEnd >> blockCount;

  if(querySize != (int) fullQuery.length()) {
    cout << "ERROR: Given PSL line doesn't fit the query sequence." << endl;
    exit(251);
  }

  DebugMsg("<ParsePSL>");
  // Parsing block size, query index and genome index strings in PSL line
  string blockSizesStr, qIndicesStr, tIndicesStr;
  str >> blockSizesStr >> qIndicesStr >> tIndicesStr;
  blockSizes = (int*) malloc (blockCount*sizeof(int));
  ParseBlockInfo(blockSizesStr, blockSizes, blockCount);
  qIndices = (int*) malloc (blockCount*sizeof(int));
  ParseBlockInfo(qIndicesStr, qIndices, blockCount);
  tIndices = (int*) malloc (blockCount*sizeof(int));
  ParseBlockInfo(tIndicesStr, tIndices, blockCount);
  DebugMsg("</ParsePSL>");

  // Gap Remove Flag is removed in basic version since restricted events require masking all gaps
  CreateAllGapRemovedQuery(fullQuery,blockCount,blockSizes,qIndices,querySize,contigDirection);

  // Accuracy calculations for the initial given alignment.
  double queryPIDinGen = calcPID(contigStart,contigEnd,genStart,genEnd,qNumInsert,matches,mismatches,repmatches);
  double alignFracQuery2Gen = calcAF(matches,mismatches,repmatches,querySize);

  for(int i=0; i<=blockCount; i++) {
    stringstream msg;
    msg << "<CheckBlock index=" << i << ">";
    DebugMsg(msg.str());
    int gapSize;
    if(i==0) {
      gapSize = qIndices[0];
    } else if(i==blockCount) {
      gapSize = fullQuery.length() - (qIndices[blockCount-1] + blockSizes[blockCount-1]);
    } else {
      gapSize = qIndices[i] - (qIndices[i-1] + blockSizes[i-1]);
    }

    if(gapSize >= gapThreshold) {
      if(gapSize < absoluteMinGapLength) {
        cout << "ERROR: there shouldn't be any gaps blatted < " << absoluteMinGapLength << "bp" << endl;
        exit(137);
      }

      int foundDuplicationFLAG = 0;

      int gapStartInd;
      if(i==0) {
        gapStartInd = 0;
      // i=blockCount case included in else
      } else {
        //for contigDirection == +
        gapStartInd = qIndices[i-1] + blockSizes[i-1];
      }

      if(contigDirection == "-") {
        //for fixing reverse listing of start indices
        gapStartInd = fullQuery.length() - gapSize - gapStartInd;
      }

      string gapSeq = fullQuery.substr(gapStartInd, gapSize);

      int numberOfNs = CountNsInGap(gapSeq);
      if(numberOfNs > (int) ((1 - queryGapMinPartialMatchThresh) * gapSize)) {
        DebugMsg("  Too many Ns\n</CheckBlock>");
        continue;
      }

      if(i==0 || i==blockCount) {
        if(CheckPolyAorT(gapSeq)) {
          DebugMsg("  PolyA\n</CheckBlock>");
          continue;
        }
      }

      // Determines the options to be used for blat calls for the gap
      string BlatOptions = globalBlatOptions;

      int gapSizeHighLimit;
      if(relaxedVersionFLAG == 0) {
        // Strict Filtering, only for duplication events (Don't look for duplications, if there aren't enough matches to align the gap)
        gapSizeHighLimit = matches + mismatches;
      } else {
        // Relaxed Filtering, sets gapsize limit allowing smallest possible duplication to be less than match
        gapSizeHighLimit =  (int) ((double) (matches + mismatches) / queryGapMinPartialMatchThresh);
      }

      if(gapSize <= gapSizeHighLimit) {
        //Create a new string with Ns in place of all gaps
        string newQuery = allGapRemovedQuery;

        char direction;
        int withinGapDupStartPos, withinGapDupEndPos, gAlignDupStartPos, gAlignDupEndPos;
        int dupMatchScore=-1;

        int eventFLAG;

        DebugMsg("<LocalAlignment>");
        if(fullContigAlignFLAG == 0) {
          eventFLAG = LocalAlignQueryWithinNeighboringArea(gapSeq, newQuery, gapStartInd, direction, withinGapDupStartPos, withinGapDupEndPos, gAlignDupStartPos, gAlignDupEndPos, dupMatchScore);
        } else {
          eventFLAG = LocalAlignQueryWithinFullContig(gapSeq, newQuery, gapStartInd, direction, withinGapDupStartPos, withinGapDupEndPos, gAlignDupStartPos, gAlignDupEndPos, dupMatchScore);
        }
        DebugMsg("</LocalAlignment>");

        if(eventFLAG == 1) {
          if(i==0 || i==blockCount) {
            //filter poly-A events on the edges even if the full gap is not poly-A
            if(CheckSubsetPolyAorT(gapSeq, withinGapDupStartPos, withinGapDupEndPos)) {
              DebugMsg("  PolyA\n</CheckBlock>");
              continue;
            }
          }

          //Fill in as similar to BLAT output
          int withinGapDupLength = withinGapDupEndPos - withinGapDupStartPos + 1;
          int gAlignDupLength = gAlignDupEndPos - gAlignDupStartPos + 1;

          foundDuplicationFLAG = 1;
          ofstream foutOut;

          //Extracting the gapSeq in forward direction even if direction is '-'
          string withinGapDupSeq = gapSeq.substr(withinGapDupStartPos, withinGapDupLength);

          //Ideal Gap Positions are introduced as original duplicated sequence indices within the contig (not within the gap)
          //these aren't modified if tandemCheck returns false
          int idealGapStartInd = gapStartInd + withinGapDupStartPos , idealGapEndInd = gapStartInd + withinGapDupEndPos;

          //Whenever ideal gap positions are refined, GAlign should be parallelly modified (or oppositely in case of inv_duplication)
          int idealGAlignStartInd = gAlignDupStartPos, idealGAlignEndInd = gAlignDupEndPos;

          int idealGapStartInd_offset=0, idealGapEndInd_offset=0;

          string ctgExtFLAG = "";
          int ctgExtReturnFLAG = 0;

          int TandemCheck = isDuplicationTandemByLocalAlign(direction, gapSize, withinGapDupLength, gAlignDupLength, gapStartInd, withinGapDupStartPos, gAlignDupStartPos, newQuery, idealGapStartInd_offset, idealGapEndInd_offset, ctgExtFLAG);

          ctgExtReturnFLAG = TandemCheck;

          idealGapStartInd += idealGapStartInd_offset;
          idealGapEndInd += idealGapEndInd_offset;
          if(direction == '+') {
            idealGAlignStartInd += idealGapStartInd_offset;
            idealGAlignEndInd += idealGapEndInd_offset;
          } else {
            idealGAlignStartInd -= idealGapEndInd_offset;
            idealGAlignEndInd -= idealGapStartInd_offset;
          }

          //Checks if non-tandem is actually caused by a tandem spanning gap (and extend within contig to restore)
          int genomicTSGFlag = 0;
          int TSGformationFlag = 0;
          if(CheckTSGFlag == 1 && direction == '+') {
            genomicTSGFlag = DoesGapOnlyAppearInContig(i, negligibleGenomicDistanceThreshold);
            if(genomicTSGFlag == 1 && TandemCheck == 0) {
              idealGapStartInd_offset=0;
              idealGapEndInd_offset=0;

              //If duplication is not touching the other edge of the gap, it is not a TSG event.
              TandemCheck = CheckTandemSpanningGap(fullQuery, gapStartInd, gapStartInd+gapSize-1, idealGapStartInd, idealGapEndInd, idealGAlignStartInd, idealGAlignEndInd, idealGapStartInd_offset, idealGapEndInd_offset);

              if(TandemCheck) {
                idealGapStartInd += idealGapStartInd_offset;
                idealGapEndInd += idealGapEndInd_offset;
                idealGAlignStartInd += idealGapStartInd_offset;
                idealGAlignEndInd += idealGapEndInd_offset;
                TSGformationFlag = 1;
              }
            }
          }

          // 1-> only Post, 2->only Pre, 3->post+pre, 4->post+novel+pre, 5->post+novel, 6->novel+pre
          int postPreFormationFlag = 0;
          int GE_postLength = 0, GE_novelLength = 0, GE_preLength =0;
          if(TandemCheck == 0 && postPreFormCheckFLAG == 1) {
            postPreFormationFlag = CheckPostPreFormation(idealGapStartInd, idealGapEndInd, idealGAlignStartInd, idealGAlignEndInd, gapStartInd, gapStartInd + gapSize - 1, fullQuery, GE_postLength, GE_novelLength, GE_preLength);
          }
          if(TandemCheck || relaxedVersionFLAG == 1) {
            foutOut.open(outputFN.c_str(),ios_base::app);
          } else {
            foutOut.open(nonBasicEventLogFile.c_str(),ios_base::app);
          }

          DebugMsg("<WriteDupEvent>");
          foutOut << "CTG:" << queryName << "(" << querySize << "bp) ";
          foutOut << "TOPOLOGY:";
          if(direction == '+') {
            if(TandemCheck) {
              foutOut << "gap-tandem-duplication ";
            } else {
              foutOut << "gap-nontandem-duplication ";
            }
          } else if(direction == '-') {
            if(TandemCheck) {
              foutOut << "gap-tandem-inverted_duplication ";
            } else {
              foutOut << "gap-nontandem-inverted_duplication ";
            }
          } else {
            cout << "ERROR: unknown direction character - Probably caused by faulty PSL line" << endl;
            exit(138);
          }

          int idealEventStartInd, idealEventEndInd, idealEventDupStartInd, idealEventDupEndInd;
          string eventMark;
          //Report event with ideal gap positions if there is inv_duplication or gap is downstream in normal duplication or duplication is non-tandem
          if(direction == '-' || idealGAlignEndInd < idealGapStartInd || TandemCheck == 0) {
            idealEventStartInd = idealGapStartInd;
            idealEventEndInd = idealGapEndInd;
            idealEventDupStartInd = idealGAlignStartInd;
            idealEventDupEndInd = idealGAlignEndInd;
            eventMark = "fromGap";
          //Report event with ideal GAlign positions if the gap is upstream in normal tandem duplication
          } else {
            idealEventStartInd = idealGAlignStartInd;
            idealEventEndInd = idealGAlignEndInd;
            idealEventDupStartInd = idealGapStartInd;
            idealEventDupEndInd = idealGapEndInd;
            eventMark = "fromGAlign";
          }

          //Checks if the inserted sequence in a non-tandem duplication is a novel seq
          int novelInsertionFLAG = 0;
          if(TandemCheck == 0) {
            int insertionStart = idealGapEndInd+1, insertionEnd = idealGAlignStartInd-1;
            if(idealGAlignEndInd < idealGapStartInd) {
              insertionStart = idealGAlignEndInd+1;
              insertionEnd = idealGapStartInd-1;
            }
            novelInsertionFLAG = CheckIfIntervalIsUnaligned(insertionStart, insertionEnd);
          }

          // Within BLAT PSL output, first base starts from 0, and for every interval (Start pos is the position of the first base, but End pos is right after the last base)
          // Furthermore the output file should consider start base as 1, and start-end positions should show inclusive interval boundaries.
          // 3rd parameters show UP end or DOWN end trimming (Alternative to N-base removal)
          int leftTrim=0, rightTrim=0;
          int idealEvent2genStart, idealEvent2genEnd;

          if(TSGformationFlag == 0) {
            idealEvent2genStart = ConvertQueryPos2GenPos(idealGAlignStartInd, leftTrim, 0);
            idealEvent2genEnd = ConvertQueryPos2GenPos(idealGAlignEndInd, rightTrim, 1);
          //There is a tandem spanning gap formation (give the genomic coordinates of the span of both duplications)
          } else {
            int firstPosToConvert = idealGapStartInd;
            if(idealGAlignStartInd < firstPosToConvert) {
              firstPosToConvert = idealGAlignStartInd;
            }

            int lastPosToConvert = idealGapEndInd;
            if(idealGAlignEndInd > lastPosToConvert) {
              lastPosToConvert = idealGAlignEndInd;
            }

            idealEvent2genStart = ConvertQueryPos2GenPos(firstPosToConvert, leftTrim, 0);
            idealEvent2genEnd = ConvertQueryPos2GenPos(lastPosToConvert, rightTrim, 1);
          }

          if(idealEvent2genStart > idealEvent2genEnd) {
            int temp = idealEvent2genStart;
            idealEvent2genStart = idealEvent2genEnd;
            idealEvent2genEnd = temp;
          }

          if(leftTrim + rightTrim > 0) {
            //Other duplication should be trimmed before reversal
            if(leftTrim>0) {
              idealEventDupStartInd += leftTrim;
            }
            if(rightTrim>0) {
              idealEventDupEndInd -= rightTrim;
            }

            if(direction == '-') {
              int temp = leftTrim;
              leftTrim = rightTrim;
              rightTrim = temp;
            }

            //gapPosition indices should be modified together with target indices
            if(leftTrim>0) {
              idealEventStartInd += leftTrim;
            }
            if(rightTrim>0) {
              idealEventEndInd -= rightTrim;
            }
          }

          int contig2genStart = genStart, contig2genEnd = genEnd-1;
          if(contigDirection == "-") {
            contig2genStart = genEnd-1;
            contig2genEnd = genStart;
          }

          foutOut << "TARGET:" << chrID << ":" << contig2genStart+1 << "-" << contig2genEnd+1 << ",";

          //if both directions are forward or both reverse, gap2gen direction is forward
          if(contigDirection.at(0) == direction) {
            foutOut << chrID << ":" << idealEvent2genStart+1 << "-" << idealEvent2genEnd+1 << " ";
          // if different, gap2gen is reverse
          } else {
            foutOut << chrID << ":" << idealEvent2genEnd+1 << "-" << idealEvent2genStart+1 << " ";
          }

          foutOut << "CONTIG:" << contigStart+1 << "-" << contigEnd << "," << idealEventStartInd+1 << "-" << idealEventEndInd+1 << " ";
          foutOut << "BREAKPOINTS:" << chrID << ":" << idealEvent2genStart+1 << "(down)-" << chrID <<":" << idealEvent2genEnd+1 << "(up)" << " ";

          double gapPIDinQuery = calcPID(withinGapDupStartPos,withinGapDupEndPos+1,gAlignDupStartPos,gAlignDupEndPos+1,0,dupMatchScore,withinGapDupLength - dupMatchScore,0);

          double alignFracGap2Query = calcAF(dupMatchScore,withinGapDupLength - dupMatchScore,0,querySize);

          char sstr[200];
          sprintf(sstr, "I1:%.2f,I2:%.2f,AF1:%.2f,AF2:%.2f", queryPIDinGen, gapPIDinQuery, alignFracQuery2Gen, alignFracGap2Query);
          foutOut << sstr << " ";

          //Block info indices for duplication events show positions within the contig
          foutOut << "BLOCKS:";
          for(int blockInd = 0; blockInd < blockCount-1; blockInd++) {
            foutOut << tIndices[blockInd]+1 << "-" << tIndices[blockInd]+blockSizes[blockInd] << ",";
          }
          foutOut << tIndices[blockCount-1]+1 << "-" << tIndices[blockCount-1]+blockSizes[blockCount-1] << ";";

          PrintGapDuplicationSplitAlignmentGenomeCoordinates(idealEvent2genStart, idealEvent2genEnd, foutOut, blockCount, blockSizes, tIndices);

          foutOut << " " << "GENES:" << genesList << " ";

          char contig2GenStrand = contigDirection.at(0);
          char gap2GenStrand = direction;
          if(contig2GenStrand == '-') {
            if(gap2GenStrand == '+') {
              gap2GenStrand = '-';
            } else {
              gap2GenStrand = '+';
            }
          }

          sprintf(sstr, "META:CS:%c,ES:%c,GAP:%d-%d,", contig2GenStrand, gap2GenStrand, gapStartInd+1, gapStartInd + gapSize);
          foutOut << sstr;

          int before_gap_pos_genome = 0, after_gap_pos_genome = 0;
          if(direction == '+') {
            before_gap_pos_genome = GetGenPos_LastBaseOfPreviousBlock(gapStartInd);
            after_gap_pos_genome = GetGenPos_FirstBaseOfFollowingBlock(gapStartInd + gapSize - 1);
          } else {
            before_gap_pos_genome = GetGenPos_LastBaseOfPreviousBlock(gapStartInd + gapSize - 1);
            after_gap_pos_genome = GetGenPos_FirstBaseOfFollowingBlock(gapStartInd);
          }

          if(contigDirection == "-") {
            int temp = before_gap_pos_genome;
            before_gap_pos_genome = after_gap_pos_genome;
            after_gap_pos_genome = temp;
          }

          if(before_gap_pos_genome == -1) {
            sprintf(sstr, "NO_GAP:N/A-N/A;%d-%d,", contig2genStart+1, contig2genEnd+1);
          } else if(after_gap_pos_genome == -1) {
            sprintf(sstr, "NO_GAP:%d-%d;N/A-N/A,", contig2genStart+1, contig2genEnd+1);
          } else {
            sprintf(sstr, "NO_GAP:%d-%d;%d-%d,", contig2genStart+1, before_gap_pos_genome+1 , after_gap_pos_genome+1 , contig2genEnd+1);
          }
          foutOut << sstr;

          int dup_Dist_Between;
          if(idealEventDupEndInd < idealEventStartInd) {
            dup_Dist_Between = idealEventStartInd - idealEventDupEndInd - 1;
          } else {
            dup_Dist_Between = idealEventDupStartInd - idealEventEndInd - 1;
          }

          sprintf(sstr, "DUP:%d-%d,DIST:%d,TRIM:%d;%d", idealEventDupStartInd+1, idealEventDupEndInd+1, dup_Dist_Between , leftTrim, rightTrim);
          foutOut << sstr;

          //Printing formation info to metadata field
          foutOut << ",FORM:";
          int numFormationsPrinted = 0;
          if(TSGformationFlag) {
            foutOut << "TSG";
            numFormationsPrinted++;
          }
          if(ctgExtReturnFLAG && ctgExtFLAG != "") {
            if(numFormationsPrinted > 0) {
              foutOut << ";";
            }
            foutOut << ctgExtFLAG;
            numFormationsPrinted++;
          }
          if(eventMark == "fromGAlign") {
            if(numFormationsPrinted > 0) {
              foutOut << ";";
            }
            foutOut << "SWAP";
            numFormationsPrinted++;
          }
          if(novelInsertionFLAG != 0) {
            if(numFormationsPrinted > 0) {
              foutOut << ";";
            }
            foutOut << "NOVEL_INS";
            if(novelInsertionFLAG == 2) {
              foutOut << "(NC)";
            }
            numFormationsPrinted++;
          }
          if(postPreFormationFlag != 0) {
            if(numFormationsPrinted > 0) {
              foutOut << ";";
            }
            // 0-> not in form, 1-> only Post, 2->only Pre, 3->post+pre, 4->post+novel+pre 5->post+novel 6->pre+novel
            switch(postPreFormationFlag)
            {
              case 1:
                foutOut << "GE_POST(" << GE_postLength << ")";
                break;
              case 2:
                foutOut << "GE_PRE(" << GE_preLength << ")";
                break;
              case 3:
                foutOut << "GE_POST_PRE(" << GE_postLength << "+" << GE_preLength << ")";
                break;
              case 4:
                foutOut << "GE_POST_NOVEL_PRE(" << GE_postLength << "+" << GE_novelLength << "+" << GE_preLength << ")";
                break;
              case 5:
                foutOut << "GE_POST_NOVEL(" << GE_postLength << "+" << GE_novelLength << ")";
                break;
              case 6:
                foutOut << "GE_NOVEL_PRE(" << GE_novelLength << "+" << GE_preLength << ")";
                break;
              default:
                cout << "ERROR: postPreFormationFLAG can only take values between 0 and 6." << endl;
                exit(137);
            }
            numFormationsPrinted++;
          }

          if(numFormationsPrinted == 0) {
            foutOut << "N/A";
          }

          if("" != extraMetaData) {
            foutOut << ",";
            foutOut << extraMetaData;
          }

          foutOut << " ";

          //When printing duplicated sequence, Ideal duplication positions are considered (not originals from BLAT alignment)
          string idealEventSequence = fullQuery.substr(idealEventStartInd, idealEventEndInd - idealEventStartInd + 1);

          foutOut << "EVENT_SEQ:" << idealEventSequence << " ";

          foutOut << "INS_SEQ:";

          if(TandemCheck) {
            foutOut << "N/A";
          } else {
            int insertIntervalStart = idealGapEndInd+1, insertIntervalEnd = idealGAlignStartInd-1;
            if(insertIntervalStart > insertIntervalEnd) {
              insertIntervalStart = idealGAlignEndInd+1;
              insertIntervalEnd = idealGapStartInd-1;
            }
            string insertSeq = fullQuery.substr(insertIntervalStart, insertIntervalEnd - insertIntervalStart + 1);
            foutOut << insertSeq;
          }
          foutOut << endl;
          DebugMsg("</WriteDupEvent>");
          numFoundMatch++;

          if(!TandemCheck) {
            numNonBasicMatch++;
          }
        }
      }

      if(foundDuplicationFLAG == 0 && gapSize >= minGapLength_InversionSearch) {
        DebugMsg("<InversionCheck>");
        // Look for inversion Events. (all for internal gaps. External gaps can be found by fusion search if long enough)
        // Firstly, in the local intronic search (only internal gaps are searched - first and last gaps aren't searched since there is no intronic boundary)
        //starting gap is i=0, ending gap is i==blockCount
        if(i>0 && i<blockCount) {
          DebugMsg("  internal block");
          int foundIntronicInversionFLAG = 0;

          int intrIndexStart = tIndices[i-1] + blockSizes[i-1];
          //end index should be +1 than actual last index than can be matched
          int intrIndexEnd = tIndices[i];

          //The intron to search should be larger
          if(intrIndexEnd - intrIndexStart >= gapSize) {
            if(tempFilesSetFLAG == 0) {
              Add_OutputDir_And_Mkstemp_To_TempFiles(outputFN);
            }
            ofstream foutGap(tempGapSeqFile.c_str());
            if(foutGap.fail()) {
              cout << "ERROR: Could not open: " << tempGapSeqFile << endl;
              exit(139);
            }

            foutGap << ">GapSeq-" << queryName << ":[" << gapStartInd << "-" << gapStartInd + gapSize - 1 << "]" << endl;
            foutGap << gapSeq << endl;
            foutGap.close();

            int BlatClientUseFlag = 0;
            if(BlatSource == "exec") {
              stringstream refWithIndices;
              refWithIndices << TwoBitGenPath << ":" << chrID << ":" << intrIndexStart << "-" << intrIndexEnd;
              string invSystemCallLine = BlatExecPath + " " + BlatOptions + " " + refWithIndices.str() + " " + tempGapSeqFile + " " + tempInvPSLOutputFile + " >> " + junkBLAT_CMDout;
              DebugMsg("<BLAT server=false>");
              system(invSystemCallLine.c_str());
              DebugMsg("</BLAT>");
            //BlatSource = "server"
            } else {
              BlatClientUseFlag = 1;
              if(BlatHostName == "LOCALHOST") {
                AssignTemporaryHostNameFile();

                stringstream systemHostCallLine;
                systemHostCallLine << "hostname >" << tempHostNameFile;
                system((systemHostCallLine.str()).c_str());
                ifstream finHostName(tempHostNameFile.c_str());
                finHostName >> BlatHostName;

                stringstream systemHostCleanLine;
                systemHostCleanLine << "rm " << tempHostNameFile;
                system((systemHostCleanLine.str()).c_str());
              }

              stringstream gfClientSystemCall;
              gfClientSystemCall << BlatClientPath << " " << BlatHostName << " " << BLAT_port_code << " " << globalBlatClientOptions << " " << BlatServer2BitDir << " " << tempGapSeqFile << " " << tempInvPSLOutputFile << " >> " << junkBLAT_CMDout;
              DebugMsg("<BLAT server=true>");
              system((gfClientSystemCall.str()).c_str());
              DebugMsg("</BLAT>");
            }

            ifstream finInvPSL(tempInvPSLOutputFile.c_str());
            int invMatch;
            while(finInvPSL >> invMatch) {
              DebugMsg("<CheckInversionMatch>");
              int invMismatch,invRepmatch,invJunk,invGapNumInsert,invGapWithinStart,invGapWithinEnd;
              finInvPSL >> invMismatch >> invRepmatch >> invJunk >> invGapNumInsert >> invJunk >> invJunk >> invJunk;
              string invDirection, invTargetName, invJunks;
              finInvPSL >> invDirection >> invJunks >> invJunk >> invGapWithinStart >> invGapWithinEnd >> invTargetName >> invJunk;

              int invPosStart, invPosEnd;
              finInvPSL >> invPosStart >> invPosEnd;

              //num blocks, block lengths, query block start indices, target block start indices.
              int invBC; string invBLs, invQSs, invTSs;
              finInvPSL >> invBC >> invBLs >> invQSs >> invTSs;

              //Discard if match isn't contained within the interval
              if(BlatClientUseFlag && (invTargetName != chrID || invPosStart < intrIndexStart || invPosEnd > intrIndexEnd)) {
                DebugMsg("  Match not in interval\n</CheckInversionMatch>");
                continue;
              }

              if(invMatch >= queryGapSimThreshold * gapSize) {
                int* invBlockSizes = (int*) malloc (invBC*sizeof(int));
                ParseBlockInfo(invBLs, invBlockSizes, invBC);
                int* invtIndices = (int*) malloc (invBC*sizeof(int));
                ParseBlockInfo(invTSs, invtIndices, invBC);

                ofstream foutOut(outputFN.c_str(),ios_base::app);

                //Blat miss cases are discarded in the Basic version
                if(invDirection == contigDirection) {
                  DebugMsg("  BLAT miss\n</CheckInversionMatch>");
                  continue;
                }

                DebugMsg("<WriteInvEvent>");
                foutOut << "CTG:" << queryName << "(" << querySize << "bp) ";
                foutOut << "TOPOLOGY:gap-internal_inversion ";

                // Gap2gen interval is assumed to be the values in the PSL line
                int gap2genStart = invPosStart;
                int gap2genEnd = invPosEnd;

                // but offset of the intronic interval is added if Blat search is done instead of gfClient
                if(BlatClientUseFlag == 0) {
                  gap2genStart += intrIndexStart;
                  gap2genEnd += intrIndexStart;

                  for(int bInd = 0; bInd<=invBC-1; bInd++) {
                    invtIndices[bInd] += intrIndexStart;
                  }
                }

                // For internal_inversion gap2gen strand is always opposite of ctg2gen
                if(contigDirection == "+") {
                  foutOut << "TARGET:" << chrID << ":" << genStart+1 << "-" << genEnd << "," << chrID << ":" << gap2genEnd << "-" << gap2genStart+1 << " ";
                } else {
                  foutOut << "TARGET:" << chrID << ":" << genEnd << "-" << genStart+1 << "," << chrID << ":" << gap2genStart+1 << "-" << gap2genEnd << " ";
                }

                int idealGapStartInd = gapStartInd + invGapWithinStart, idealGapEndInd= gapStartInd + invGapWithinEnd - 1;

                foutOut << "CONTIG:" << contigStart+1 << "-" << contigEnd << "," << idealGapStartInd+1 << "-" << idealGapEndInd+1 << " ";
                foutOut << "BREAKPOINTS:" << chrID << ":" << gap2genStart+1 << "(down)-" << chrID <<":" << gap2genEnd << "(up)" << " ";

                double invGapPIDinQuery = calcPID(invGapWithinStart,invGapWithinEnd,invPosStart,invPosEnd,invGapNumInsert,invMatch,invMismatch,invRepmatch);
                double invAlignFracGap2Query = calcAF(invMatch,invMismatch,invRepmatch,querySize);
                char sstr[200];
                sprintf(sstr, "I1:%.2f,I2:%.2f,AF1:%.2f,AF2:%.2f", queryPIDinGen, invGapPIDinQuery, alignFracQuery2Gen, invAlignFracGap2Query);
                foutOut << sstr << " ";

                //Block info for intron-wide inversion events show positions within cropped genome (intron within neighbour blocks)
                foutOut << "BLOCKS:";
                for(int blockInd = 0; blockInd < blockCount-1; blockInd++) {
                  foutOut << tIndices[blockInd]+1 << "-" << tIndices[blockInd]+blockSizes[blockInd] << ",";
                }
                foutOut << tIndices[blockCount-1]+1 << "-" << tIndices[blockCount-1]+blockSizes[blockCount-1] << ";";

                //Since target block info is kept in forward direction in BLAT PSL format, direction doesn`t matter
                for(int blockInd = 0; blockInd < invBC-1; blockInd++) {
                  int g2genomeIndexStart = invtIndices[blockInd];
                  int g2genomeIndexEnd = invtIndices[blockInd] + invBlockSizes[blockInd] - 1;
                  foutOut << g2genomeIndexStart+1 << "-" << g2genomeIndexEnd+1 <<  ",";
                }
                int g2genomeIndexStart = invtIndices[invBC-1];
                int g2genomeIndexEnd = invtIndices[invBC-1] + invBlockSizes[invBC-1] - 1;
                foutOut  << g2genomeIndexStart+1 << "-" << g2genomeIndexEnd+1 << " ";

                foutOut << "GENES:" << genesList << " ";

                char contig2GenStrand = contigDirection.at(0), gap2GenStrand = invDirection.at(0);

                if(contigDirection == "+") {
                  sprintf(sstr,"META:CS:%c,ES:%c,GAP:%d-%d,NO_GAP:%d-%d;%d-%d,DUP:N/A-N/A,DIST:N/A,TRIM:N/A;N/A", contig2GenStrand, gap2GenStrand, gapStartInd+1, gapStartInd + gapSize, genStart+1, intrIndexStart - 1, intrIndexEnd, genEnd);
                } else {
                  sprintf(sstr,"META:CS:%c,ES:%c,GAP:%d-%d,NO_GAP:%d-%d;%d-%d,DUP:N/A-N/A,DIST:N/A,TRIM:N/A;N/A", contig2GenStrand, gap2GenStrand, gapStartInd+1, gapStartInd + gapSize, genEnd, intrIndexEnd, intrIndexStart - 1, genStart+1);
                }
                foutOut << sstr;

                foutOut << ",FORM:N/A";

                if("" != extraMetaData) {
                  foutOut << ",";
                  foutOut << extraMetaData;
                }

                foutOut << " ";

                // Ideal gap relocation isn't considered for internal inversions (it requires extension through genome).
                // Different than duplication EVENT_SEQ, here it represents the sequence starting from alignment start and end points within the original gap interval;
                // and doesn't give any info related to multiple blocks.
                string eventSequenceSimple = gapSeq.substr(invGapWithinStart, invGapWithinEnd - invGapWithinStart);
                foutOut << "EVENT_SEQ:" << eventSequenceSimple << " ";
                foutOut << "INS_SEQ:N/A" << endl;

                free(invBlockSizes);
                free(invtIndices);

                DebugMsg("</WriteInvEvent>");
                numFoundMatch++;
                foundIntronicInversionFLAG = 1;
              }
              DebugMsg("</CheckInversionMatch>");
            }
            finInvPSL.close();
          }
        } else {
          DebugMsg("  edge block");
        }
        DebugMsg("</InversionCheck>");
      }
    }
    DebugMsg("</CheckBlock>");
  }
  DebugMsg("</GapSearch>");
  return numFoundMatch;
}

void CheckAndUpdateConfig(string configFileName) //Reads a config file and replaces with the default path and options if they exist.
{
  ifstream fin(configFileName.c_str());
  if(fin.is_open()) {
    string line;
    while(getline(fin,line))
    {
      if(line.length()>11 && line.substr(0,11) == "BLATsource=") {
        BlatSource = line.substr(11, line.length()-11);
      }
      if(line.length()>9 && line.substr(0,9) == "BLATpath=") {
        BlatExecPath = line.substr(9,line.length()-9);
      }
      if(line.length()>13 && line.substr(0,13) == "gfClientPath=") {
        BlatClientPath = line.substr(13,line.length()-13);
      }
      if(line.length()>16 && line.substr(0,16) == "BLATserver_host=") {
        BlatHostName = line.substr(16,line.length()-16);
      }
      if(line.length()>19 && line.substr(0,19) == "BLATserver_2bitDir=") {
        BlatServer2BitDir = line.substr(19,line.length()-19);
      }
      if(line.length()>12 && line.substr(0,12) == "BLAToptions=") {
        globalBlatOptions = line.substr(12,line.length()-12);
      }
      if(line.length()>16 && line.substr(0,16) == "gfClientOptions=") {
        globalBlatClientOptions = line.substr(16,line.length()-16);
      }
      if(line.length()>12 && line.substr(0,12) == "2bitGenPath=") {
        TwoBitGenPath = line.substr(12,line.length()-12);
      }
      if(line.length()>22 && line.substr(0,22) == "2bitToFaConverterPath=") {
        TwoBitToFaConverterPath = line.substr(22, line.length()-22);
      }
      if(line.length()>15 && line.substr(0,15) == "relaxedVersion=") {
        if(line.at(15) == '1') {
          relaxedVersionFLAG = 1;
        }
      }
      if(line.length()>16 && line.substr(0,16) == "fullContigAlign=") {
        if(line.at(16) == '1') {
          fullContigAlignFLAG = 1;
        }
      }
      if(line.length()>19 && line.substr(0,19) == "dupAlignGapPenalty=") {
        optionalGapPenalty = atoi((line.substr(19,line.length()-19)).c_str());
      }
      if(line.length()>9 && line.substr(0,9) == "TSGcheck=") {
        CheckTSGFlag = atoi((line.substr(9,line.length()-9)).c_str());
      }
      if(line.length()>10 && line.substr(0,10) == "genomeExt=") {
        postPreFormCheckFLAG = atoi((line.substr(10,line.length()-10)).c_str());
      }
    }
  } else {
    cout << "ERROR: could not open gap config file: \"" << configFileName << "\"." << endl;
    exit(245);
    }

  //Error checking
  if(BlatSource == "server") {
    if(BlatHostName == "" || BlatServer2BitDir == "" || globalBlatClientOptions == "" || BlatClientPath == "") {
      cout << "ERROR: gfServer/gfClient use for inversion search requires BLATserver_host, BLATserver_2bitDir, gfClientPath, and gfClientOptions with required -nohead option" << endl;
      exit(242);
    }
  } else if(BlatSource == "exec") {
    if(BlatExecPath == "" || TwoBitGenPath == "" || globalBlatOptions == "") {
      cout << "ERROR: BLATsource=exec option requires both BLAT executable path \"BLATpath=\" and 2bit genome file path \"2bitGenPath=\" and BlatOptions with required -noHead option" << endl;
      exit(243);
    }
  } else {
    cout << "ERROR: Blat Source must be \"server\" or \"exec\". \"" << BlatSource << "\" is not a valid Blat Source." << endl;
    exit(241);
  }

  if(postPreFormCheckFLAG == 1) {
    if(TwoBitGenPath == "" || TwoBitToFaConverterPath == "") {
      cout << "ERROR: Detecting genomic extension forms requires both 2bit genome path \"2bitGenPath=\" and 2bit to fasta converter path \"2bitToFaConverterPath=\"." << endl;
      exit(244);
    }
  }
}

void PrintProgramSpecifications() {
  cout << "================================================================================\n"
       << "Gap Realigner " << ver <<" specifications are as follows:                       \n"
       << "================================================================================\n"
       << "Call with '-v' for version, '-F' for formation output field types,              \n"
       << " -E for exit codes                                                              \n"
       << "================================================================================\n"
       << "Required Command Line Arguments:                                                \n"
       << "1) [PSL_LINE]         : A full PSL line within quotation marks                  \n"
       << "2) [CONTIG_SEQ]       : The original sequence that PSL line is obtained from    \n"
       << "3) [GENE_LIST_STRING] : A string of genes to be reported in the output          \n"
       << "4) [OUTPUT_FILENAME]  : Name of the output file (reported lines are appended)   \n"
       << "5) [CONFIG_FILENAME]  : Name of the configuration file (specifications below)   \n"
       << "--------------------------------------------------------------------------------\n"
       << "Optional Parameters:                                                            \n"
       << "6) [Gap Length Threshold]: [>=3] Minimum length to consider as a gap            \n"
       << "7) [Similarity Threshold]: [ 0.67 <= double <= 1.0]                             \n"
       << "      Minimum percentage identity threshold for aligned regions                 \n"
       << "8) [Partial Alignment Ratio]: [ 0 < double <= 1]                                \n"
       << "      Minimum ratio of partial alignment to the overall gap                     \n"
       << "9) [Extra Meta Data]: A string that will be reported at the end of the meta data\n"
       << "10) [BLAT Server Port No]: Port number for gfClient                             \n"
       << "      (required if BLATsource=server)                                           \n"
       << "================================================================================\n"
       << "Configuration File Specifications:                                              \n"
       << "--------------------------------------------------------------------------------\n"
       << "   *)  Parameters can be written in any order                                   \n"
       << "  **)  Identifiers should start from beginning of the line and should not       \n"
       << "       contain any blankspace other than option parameters                      \n"
       << "       (BLAToptions and gfClientOptions)                                        \n"
       << " ***)  Should be written with '=' between the identifier and parameter          \n"
       << "       (e.g. BLATsource=exec or 2bitGenPath=hg18.2bit)                          \n"
       << "--------------------------------------------------------------------------------\n"
       << "BLATsource            : [exec/server] Whether BLAT is an executable or gfServer \n"
       << "BLATpath              : Path for BLAT executable                                \n"
       << "BLAToptions           : List of options for BLAT calls (\"-noHead\" required)   \n"
       << "2bitGenPath           : Path for 2bit genome file                               \n"
       << "2bitToFaConverterPath : Path for TwoBitToFa executable (required if genomeExt=1)\n"
       << "relaxedVersion        : [0/1 : Default(0)]                                      \n"
       << "  0 ensures strict match ratio similarity                                       \n"
       << "  1 allows relaxed partial mapping similarity thresholds                        \n"
       << "fullContigAlign       : [0/1 : Default(0)] when looking for duplications        \n"
       << "                        if 1, gap is aligned to full sequence                   \n"
       << "                        if 0, the gap is aligned to only within 2x vicinity     \n"
       << "BLATserver_host       : Hostname for gfServer                                   \n"
       << "  (LOCALHOST for local machine adress)                                          \n"
       << "BLATserver_2bitDir    : Directory of the 2bit file used by BLAT gfServer        \n"
       << "                        (required if BLATsource=server)                         \n"
       << "gfClientPath          : path for gfClient executable                            \n"
       << "  (required if BLATsource=server)                                               \n"
       << "gfClientOptions       : options for gfClient executable (\"-nohead\" required)  \n"
       << "                        (required if BLATsource=server)                         \n"
       << "dupAlignGapPenalty    : [positive int : Default(1)]                             \n"
       << "                        Constant gap penalty value for duplication SW-alignment \n"
       << "TSGcheck              : [0/1 : Default(1)]                                      \n"
       << "  Whether tandem spanning gaps are checked                                      \n"
       << "genomeExt             : [0/1 : Default(0)]                                      \n"
       << "  Determines if Post-Novel-Pre insertion forms are to be detected               \n"
       << "================================================================================\n"
       << endl;
}

void CheckArgumentValues(int portCode)
{
  if(BlatSource == "server" && (portCode < 100 || portCode > 70000)) {
    PrintProgramSpecifications();
    cout << "ERROR: Specified port code " << portCode << " isn't within port range" << endl;
    exit(252);
  }
  if(gapThreshold < absoluteMinGapLength) {
    PrintProgramSpecifications();
    cout << "ERROR:   " <<  absoluteMinGapLength << "<= Gap Length Threshold" << endl;
    exit(253);
  }
  if(queryGapSimThreshold < 0.67 || queryGapSimThreshold > 1.0) {
    PrintProgramSpecifications();
    cout << "ERROR: 0.67 <= Similarity Threshold <= 1.0" << endl;
    exit(254);
  }
  if(queryGapMinPartialMatchThresh <=0 || queryGapMinPartialMatchThresh > 1.0) {
    PrintProgramSpecifications();
    cout << "ERROR: 0 < Partial Alignment Ratio <= 1" << endl;
    exit(255);
  }
}

void PrintFormationTypes()
{
  cout << "================================================================================\n"
       << "  Gap Realigner " << ver << " formation types are as follows:                   \n"
       << "================================================================================\n"
       << "TSG: Tandem spanning gap.                                                       \n"
       << "  Gap originally appears as spanning a tandem forward duplication.              \n"
       << "  Event is reported as the contig downstream copy of the duplication.           \n"
       << "--------------------------------------------------------------------------------\n"
       << "SWAP: Ideal event is swapped                                                    \n"
       << "  If the original gap in tandem forward duplication is on the upstream side,    \n"
       << "  event coordinates are swapped with downstream copy.                           \n"
       << "--------------------------------------------------------------------------------\n"
       << "TAN([type]): Reports how tandemness is detected, [type] can be:                 \n"
       << "  Sh&Neg : The distance between the duplicated regions is <2bp                  \n"
       << "           Borders have been modified accordingly.                              \n"
       << "           This is the formation reported for an ideal tandem duplication.      \n"
       << "  Sh&Ir  : The distance between the duplicated regions was short in length      \n"
       << "           relative to the duplications. However a regular extension would      \n"
       << "           exceed the contig border, thus returned as tandem but borders        \n"
       << "           aren't modified.                                                     \n"
       << "  Sh&NMod: The distance between duplicated regions is short enough to warrant   \n"
       << "           for tandemness, but borders aren't modified since the extending      \n"
       << "           sequences were not similar enough.                                   \n"
       << "  Sh&Mod : Same as above; but this time extending sequences were similar also,  \n"
       << "           thus borders of the duplication are modified.                        \n"
       << "  NMod   : Distance between the duplicated regions were longer, but the         \n"
       << "           sequences could be extended to make them tandem, and the overall     \n"
       << "           extended sequences were similar enough to accept as tandem           \n"
       << "           duplication. However extending sequences themselves weren't similar  \n"
       << "           enough, thus the borders are not modified.                           \n"
       << "  Mod    : Same as above, but in this case the extending sequences are similar  \n"
       << "           also thus borders of the duplication are modified.                   \n"
       << "--------------------------------------------------------------------------------\n"
       << "NOVEL_INS: Novel insertion form                                                 \n"
       << "       In a non-tandem duplication event, indicates that the insertion between  \n"
       << "       the duplicated sequences has not been aligned to the genome.             \n"
       << "       If it appears as NOVEL_INS(NC), indicates that majority of insertion is  \n"
       << "       unaligned yet there is a negligibly short overlap on either side with    \n"
       << "       another block (This can happen due to similarity of several bases near   \n"
       << "       the border) This overlap length should be at most:                       \n"
       << "           1 + (int) ( insertion_length * (1 - queryGapSimThreshold) / 2.0)     \n"
       << "--------------------------------------------------------------------------------\n"
       << "GE_[type]([lengths]): Genomic extension forms in non-tandem duplications        \n"
       << "       POST          : The insertion between two duplications corresponds to    \n"
       << "                       the genomic region right after the duplicated sequence.  \n"
       << "       PRE           : The insertion between the duplications corresponds to    \n"
       << "                       the genomic region right before the duplicated sequence. \n"
       << "       POST_NOVEL    : Same as GE_POST formation but in this case there is a    \n"
       << "                       novel insertion in the contig between the POST region    \n"
       << "                       and downstream duplication.                              \n"
       << "       NOVEL_PRE     : Same as GE_PRE but there is a novel insertion in contig  \n"
       << "                       between upstream duplication and PRE region.             \n"
       << "       POST_PRE      : The insertion between two duplications corresponds to    \n"
       << "                       the combination of sequences from POST and PRE regions   \n"
       << "                       in the genome.                                           \n"
       << "       POST_NOVEL_PRE: Same as GE_POST_PRE but there is a novel insertion in    \n"
       << "                       the contig between the POST and PRE regions.             \n"
       << "                                                                                \n"
       << "       For all of these types the length of these genomic extensions or novel   \n"
       << "       sequences are given right after the formation type name                  \n"
       << "       (e.g.: POST(9)    NOVEL_PRE(4+22)     POST_NOVEL_PRE(10+5+12) )          \n"
       << "================================================================================\n"
       << endl;
}

void PrintExitCodes()
{
  cout << "================================================================================\n"
       << "  Gap Realigner " << ver << " exit codes are as follows:                        \n"
       << "================================================================================\n"
       << " Regular exit codes:                                                            \n"
       << "--------------------------------------------------------------------------------\n"
       << " [0]   : No gap events are detected                                             \n"
       << " [1-99]: Number of detected gap events is the returned number                   \n"
       << " [100] : Number of detected gap events is >= 100                                \n"
       << " [110] : In the un-relaxed version,                                             \n"
       << "         all of the caught gap events are non-basic events                      \n"
       << " [250] : Program exits after printing program specifications, version info,     \n"
       << "         formation description, or exit code information.                       \n"
       << "================================================================================\n"
       << " Input error exit codes:                                                        \n"
       << "--------------------------------------------------------------------------------\n"
       << " [251-5] : Argument Parameter Faults (one of the following)                     \n"
       << "   [251] - Given PSL line doesn't fit the query sequence                        \n"
       << "   [252] - Specified port code isn't in range                                   \n"
       << "   [253] - Gap length threshold is not non-negative                             \n"
       << "   [254] - Similarity Threshold is larger than 1 or smaller than 0.67           \n"
       << "   [255] - Partial alignment ration is not non-negative or larger than 1        \n"
       << "--------------------------------------------------------------------------------\n"
       << " [241-5] : Configuration parameter faults (one of the following)                \n"
       << "   [241] - BLAT source isn't determined as executable (exec) or                 \n"
       << "           server/client (server)                                               \n"
       << "   [242] - BLAT source option is selected as executable, but BLAT executable    \n"
       << "           path or 2bit genome file path aren't determined. Or BlatOptions line \n"
       << "           doesn't contain the required -noHead option.                         \n"
       << "   [243] - BLAT source option is selected as server, however hostname, client   \n"
       << "           executable path or BLATserver_2bitDir aren't given as described in   \n"
       << "           the specifications. Or gfClientOptions line doesn't contain the      \n"
       << "           required -nohead option. (Note that -nohead is for client options,   \n"
       << "           -noHead is for Blat options.)                                        \n"
       << "   [244] - If Post-Pre formation checking is on, 2bit genome path and           \n"
       << "           2bitToFasta converter executable path aren't provided.               \n"
       << "   [245] - Could not open gap config file.                                      \n"
       << "================================================================================\n"
       << " Runtime error exit codes:                                                      \n"
       << "--------------------------------------------------------------------------------\n"
       << " [135] : There is an unexpected character in input sequence, or extracted       \n"
       << "         genomic sequence. Acceptable characters are {A,C,G,T,N,a,c,g,t,n}      \n"
       << " [136] : Length of the compared sequences in OneSidedLocalAlignment function    \n"
       << "         are not equal. This indicates an error in Post-Pre formation functions.\n"
       << " [137] : A full gap to be aligned has less than 7 bases. Indicates an error     \n"
       << "         with gap length filtering or similarity/partial ratio calculation      \n"
       << " [138] : The direction character is not {+,-}                                   \n"
       << "         This might be caused by a faulty PSL line input.                       \n"
       << " [139] : Temporary gap sequence file could not be opened. Might be caused by a  \n"
       << "         diskspace issue or Mkstemp function access                             \n"
       << " [140] : Ideal gap and ideal gap alignment (gAlign) intervals are overlapping,  \n"
       << "         where as they should be seperate intervals. This might indicate an     \n"
       << "         error in Tandem spanning gap detection functions.                      \n"
       << " [141] : While converting a give query position to genomic position. The query  \n"
       << "         index is larger than the query size. This indicates an query pos       \n"
       << "         calculation error.                                                     \n"
       << "================================================================================\n"
       << endl;
}

int main(int argc, char *argv[])
{
  if(1 < argc && string("-v") == string(argv[1])) {
    cout << "Gap Realigner Version: " << ver << endl;
    exit(250);
  }

  if(1 < argc && string("-F") == string(argv[1])) {
    PrintFormationTypes();
    exit(250);
  }

  if(1 < argc && string("-E") == string(argv[1])) {
    PrintExitCodes();
    exit(250);
  }

  if(argc < 6) {
    PrintProgramSpecifications();
    exit(250);
  }

  // PSL line
  stringstream PSL_str(argv[1]);
  // Full query sequence
  string fullQuery(argv[2]);
  // List of genes to report (not used in gap realigner)
  string gene_list(argv[3]);
  // Path of the output file
  string outputFileName(argv[4]);
  // Path of the configration file
  string blatConfigFile(argv[5]);

  int blatServerPortCode = 0;
  if(argc > 6) {
    //Length threshold for a gap to be considered
    gapThreshold = atoi(argv[6]);
  }
  if(argc > 7) {
    //Similarity threshold for a gap or partial gap alignment to be considered as successful
    queryGapSimThreshold = atof(argv[7]);
  }
  if(argc > 8) {
    //Ratio threshold of a partial gap to the overall gap to be considered as valid
    queryGapMinPartialMatchThresh = atof(argv[8]);
  }
  if(argc > 9) {
    //Extra string to include at end of meta data in output
    extraMetaData = argv[9];
  }
  if(argc > 10) {
    if (string("-d") == string(argv[10])) {
      debug = true;
      if(argc > 11) {
        //Port number to be used in BLAT server/client interaction
        blatServerPortCode = atoi(argv[11]);
      }
    } else {
      //Port number to be used in BLAT server/client interaction
      blatServerPortCode = atoi(argv[11]);
    }
  }

  DebugMsg("<Main>");
  //Loads configuration parameter values from the config file
  CheckAndUpdateConfig(blatConfigFile);
  // Checks if argument values and config parameter values are valid
  CheckArgumentValues(blatServerPortCode);
  if(postPreFormCheckFLAG == 1) {
    // Determines random file extensions for temporary files
    Add_OutputDir_And_Mkstemp_To_TempFiles(outputFileName);
  }

  //Search for a gap event
  int numFound = GapSearch(PSL_str, fullQuery, outputFileName, gene_list, blatServerPortCode);

  if(tempFilesSetFLAG) {
    CleanTempFiles();
  }

  //If all found events are non-basic events in the un-relaxed version of the code, don't return number, return errorCode instead
  if(relaxedVersionFLAG == 0 && 0 < numFound && numFound == numNonBasicMatch) {
    DebugMsg("<OnlyNonBasic/>\n</Main>");
    exit(110);
  }

  //Doesn't return anything larger than 100, to leave the rest for error exits.
  if(numFound > 100) {
    DebugMsg("<MoreThan100/>\n</Main>");
    numFound = 100;
  }

  //Returns the number of detected gap events as output
  DebugMsg("<Normal/>\n</Main>");
  exit(numFound);

  return 0;
}

