import os
import logging
import subprocess
logging.basicConfig(level=logging.INFO)


class CondorHandler(object) :
  """
    A class to submit batch jobs to a HTCondor scheduler.

    ...

    Attributes
    ----------
    batch_path : str
        path where the batch config files which are created will be stored
    log_path : str
        path where the batch log files will be stored


    Methods
    -------
    activate_testmode():
        Activate test mode, allowing for checks of config files in dry runs, no jobs submitted
    deactivate_testmode():
        Dctivate test mode, enable submitting jobs
    send_job(command, tag = "htcondor_job"):
        Submit job by creating config files (bash file and HTCondor submission file) and executing condor_submit

  """
  def __init__(self, batch_path, log_path):
    self.batch_path = batch_path
    self.log_path = log_path

    # internal options
    self._tag = "htcondor_job"
    self._condor_options = {}
    self._test_mode = False


  def __setitem__(self, key, value):
    self._condor_options[key] = value


  def __getitem__(self, key):
    return self._condor_options[key]


  def activate_testmode(self):
    logging.debug("Activated test mode: not submitting any jobs.")
    self._test_mode = True


  def deactivate_testmode(self):
    logging.debug("Deactivated test mode: submitting jobs.")
    self._test_mode = False


  def send_job(self, command, tag="htcondor_job"):
    # make files
    self._tag = tag
    bashfile = self._make_bash_file(command)
    jobfile = self._make_job_file(bashfile)
    # do submit thing
    if self._test_mode:
      logging.debug("Created job file {0}".format(jobfile))
    else:
      subprocess.call("condor_submit {0}".format(jobfile), shell=True)


  def _make_bash_file(self, command) :
    runFile = os.path.join(self.batch_path, "batch_{0}.sh".format(self._tag))
    with open(runFile,"w") as fr :
      fr.write('#!/bin/sh\n')
      fr.write('# {0} batch run script\n'.format(self._tag))
      fr.write('#$ -cwd\n')
      fr.write('#$ -j y\n')
      fr.write('#$ -l cvmfs\n')
      fr.write("path=" + os.getcwd() + "\n")
      fr.write('export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase\n')
      fr.write('source /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/user/atlasLocalSetup.sh\n')
      fr.write('pwd; ls -l\n')
      fr.write(command + "\n")
      fr.write('ls -l\n')
      fr.close()
    os.system("chmod a+x " + runFile)
    logging.debug("Made run file " + runFile)
    return runFile


  def _make_job_file(self, runFile) :
    batchFile = os.path.join(self.batch_path, "batch_{0}.job".format(self._tag))
    with open(batchFile, "w") as fs :
      if 'universe' in self._condor_options.keys():
        fs.write('Universe            = {0}\n'.format(self._condor_options['universe']))
      if 'jobflavour' in self._condor_options.keys():
        fs.write('+JobFlavour         = "{0}"\n'.format(self._condor_options['jobflavour']))
      if 'project' in self._condor_options.keys():
        fs.write("+MyProject          = \"{0}\"\n".format(self._condor_options['project']))
      if 'runtime' in self._condor_options.keys():
        fs.write("+RequestRuntime     = {0}\n".format(self._condor_options['runtime']))
      if 'memory' in self._condor_options.keys():
        fs.write("Request_Memory      = {0}\n".format(self._condor_options['memory']))
      if 'cpu' in self._condor_options.keys():
        fs.write("Request_CPUs        = {0}\n".format(self._condor_options['cpu']))
      if 'requirements' in self._condor_options.keys():
        fs.write("Requirements        = {0}\n".format(self._condor_options['requirements']))
      if 'container' in self._condor_options.keys():
        fs.write("+MySingularityImage = \"{0}\"\n".format(self._condor_options['container']))
      fs.write('Executable          = {0}\n'.format(runFile))
      fs.write('Output              = {0}/stdout_{1}_$(ClusterId).txt\n'.format(self.log_path, self._tag))
      fs.write('Error               = {0}/stderr_{1}_$(ClusterId).txt\n'.format(self.log_path, self._tag))
      fs.write('log                 = {0}/batch_{1}_$(ClusterId).log\n'.format(self.log_path, self._tag))
      fs.write('\nqueue\n')
      fs.close()
    logging.debug("Made job file " + batchFile)
    return batchFile

