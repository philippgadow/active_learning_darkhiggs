/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

/***********************************************************************/
/*                                                                     */
/*              RestFramesHelper                                       */
/*                                                                     */
/*              Author: Christopher Rogan                              */
/*              April, 2017                                            */
/*                                                                     */
/***********************************************************************/
#include "SimpleAnalysisFramework/RestFramesHelper.h"

RestFramesHelper::~RestFramesHelper() {
  int N = m_Objects.size();
  for (int i = 0; i < N; i++)
    delete m_Objects[i];

  m_Objects.clear();
  m_LabFrames.clear();
  m_DecayFrames.clear();
  m_SAFrames.clear();
  m_VisFrames.clear();
  m_InvFrames.clear();
  m_InvJigsaws.clear();
  m_CombJigsaws.clear();
  m_InvGroups.clear();
  m_CombGroups.clear();
}

LabRecoFrame *RestFramesHelper::addLabFrame(const std::string &name) {
  LabRecoFrame *frame = new LabRecoFrame(name, name);

  m_Objects.push_back(frame);
  m_LabFrames[name] = frame;

  return frame;
}

DecayRecoFrame *RestFramesHelper::addDecayFrame(const std::string &name) {
  DecayRecoFrame *frame = new DecayRecoFrame(name, name);

  m_Objects.push_back(frame);
  m_DecayFrames[name] = frame;

  return frame;
}

VisibleRecoFrame *RestFramesHelper::addVisibleFrame(const std::string &name) {
  VisibleRecoFrame *frame = new VisibleRecoFrame(name, name);

  m_Objects.push_back(frame);
  m_VisFrames[name] = frame;

  return frame;
}

InvisibleRecoFrame *
RestFramesHelper::addInvisibleFrame(const std::string &name) {
  InvisibleRecoFrame *frame = new InvisibleRecoFrame(name, name);

  m_Objects.push_back(frame);
  m_InvFrames[name] = frame;

  return frame;
}

SelfAssemblingRecoFrame *RestFramesHelper::addSAFrame(const std::string &name) {
  SelfAssemblingRecoFrame *frame = new SelfAssemblingRecoFrame(name, name);

  m_Objects.push_back(frame);
  m_SAFrames[name] = frame;

  return frame;
}

InvisibleJigsaw *RestFramesHelper::addInvisibleJigsaw(const std::string &name,
                                                      InvJigsawType type) {
  InvisibleJigsaw *jigsaw = nullptr;

  if (type == kSetMass)
    jigsaw = new SetMassInvJigsaw(name, name);

  if (type == kSetRapidity)
    jigsaw = new SetRapidityInvJigsaw(name, name);

  if (type == kContraBoost)
    jigsaw = new ContraBoostInvJigsaw(name, name);

  m_Objects.push_back(jigsaw);
  m_InvJigsaws[name] = jigsaw;

  return jigsaw;
}

MinMassesCombJigsaw *
RestFramesHelper::addCombinatoricJigsaw(const std::string &name,
                                        CombJigsawType type) {
  MinMassesCombJigsaw *jigsaw = nullptr;

  if (type == kMinMasses)
    jigsaw = new MinMassesCombJigsaw(name, name);

  m_Objects.push_back(jigsaw);
  m_CombJigsaws[name] = jigsaw;

  return jigsaw;
}

InvisibleGroup *RestFramesHelper::addInvisibleGroup(const std::string &name) {
  InvisibleGroup *group = new InvisibleGroup(name, name);

  m_Objects.push_back(group);
  m_InvGroups[name] = group;

  return group;
}

CombinatoricGroup *
RestFramesHelper::addCombinatoricGroup(const std::string &name) {
  CombinatoricGroup *group = new CombinatoricGroup(name, name);

  m_Objects.push_back(group);
  m_CombGroups[name] = group;
  return group;

  return group;
}

LabRecoFrame *RestFramesHelper::getLabFrame(const std::string &name) {
  return m_LabFrames[name];
}

DecayRecoFrame *RestFramesHelper::getDecayFrame(const std::string &name) {
  return m_DecayFrames[name];
}

VisibleRecoFrame *RestFramesHelper::getVisibleFrame(const std::string &name) {
  return m_VisFrames[name];
}

InvisibleRecoFrame *
RestFramesHelper::getInvisibleFrame(const std::string &name) {
  return m_InvFrames[name];
}

SelfAssemblingRecoFrame *RestFramesHelper::getSAFrame(const std::string &name) {
  return m_SAFrames[name];
}

InvisibleJigsaw *RestFramesHelper::getInvisibleJigsaw(const std::string &name) {
  return m_InvJigsaws[name];
}

MinMassesCombJigsaw *
RestFramesHelper::getCombinatoricJigsaw(const std::string &name) {
  return m_CombJigsaws[name];
}

InvisibleGroup *RestFramesHelper::getInvisibleGroup(const std::string &name) {
  return m_InvGroups[name];
}

CombinatoricGroup *
RestFramesHelper::getCombinatoricGroup(const std::string &name) {
  return m_CombGroups[name];
}
