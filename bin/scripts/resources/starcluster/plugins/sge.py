#!/usr/bin/env python
import os
import re
import string
import sys
import time
import posixpath
import subprocess


from starcluster.clustersetup import ClusterSetup
from starcluster.logger import log
from starcluster import utils

class NullDevice():
    def write(self, s):
        pass

class CreateCell (ClusterSetup):

    """
    Configures a custom SGE cell for a StarCluster cluster
    """
    def __init__(self, head_ip, head_long_ip, cell, execdport, qmasterport, root, slots):
        log.info("Loaded plugin: sge.CreateCell")

        log.debug("CreateCell.__init__    Initialising CreateCell plugin.")
        log.debug("CreateCell.__init__    head_ip %s" % head_ip)        
        log.debug("CreateCell.__init__    head_long_ip %s" % head_long_ip)        
        log.debug("CreateCell.__init__    cell %s" % cell)
        log.debug("CreateCell.__init__    execdport %s" % execdport)
        log.debug("CreateCell.__init__    qmasterport %s" % qmasterport)
        log.debug("CreateCell.__init__    root %s" % root)
        log.debug("CreateCell.__init__    slots %s" % slots)

        self.head_ip = head_ip
        self.head_long_ip = head_long_ip
        self.cell = cell
        self.execdport = execdport
        self.qmasterport = qmasterport
        self.root = root

        os.environ['SGE_ROOT'] = root
        os.environ['SGE_CELL'] = cell
        os.environ['SGE_QMASTER_PORT'] = qmasterport
        os.environ['SGE_EXECD_PORT'] = execdport


    def run(self, nodes, master, user, user_shell, volumes):
        """
            Mount NFS shares on master and all nodes
        """
        log.info("Running plugin: sge.CreateCell")
        log.debug("CreateCell.run    CreateCell.run(nodes, master, user, user_shell, volumes")

        #### OPEN NEW PORTS ON EC2 ON HEAD
        self.openSgePorts()

        #### ADD ENVIRONMENT VARIABLES TO /etc/profile ON MASTER/NODES
        for node in nodes:
           self.addEnvarsToProfile(node)

        #### CREATE NEW CELL DIRECTORY ON HEAD AND MASTER/NODES
        self.copyCellOnHead()
        #for node in nodes:
        self.copyCell(master)

        #### UPDATE HEAD'S IP ADDRESS (MAY HAVE CHANGED IF REBOOTED)
        self.updateHeadIp()

        #### UPDATE MASTER'S IP ADDRESS IN act_qmaster FILE ON HEAD
        self.updateMasterIp(master);

        ##### RESTART SGE ON MASTER/NODES
        for node in nodes:
            self.restartSge(node)

        ####  SCHEDULING INFO
        self.enableSchedulingInfo()

        #### ADD threaENABLEded PARALLEL ENVIRONMENT ON MASTER
        self.addParallelEnvironment(master)

        #### ADD HEAD TO SUBMIT AND ADMIN HOSTS ON MASTER
        self.addHeadToHosts(master)

        #### ADD NODES TO @allhosts GROUP
        for node in nodes:
            self.addToAllhosts(node, master)

        ##### RESTART SGE ON MASTER/NODES
        for node in nodes:
            self.restartSge(node)

        log.info("Completed running plugin: sge.CreateCell")

    def addToAllhosts(self, node, master):
        """
            Add host to @allhosts group to enable it to be an execution host
        """
        log.info("Add %s to @allhosts group", node.alias)

        os.environ['SGE_ROOT'] = self.root
        os.environ['SGE_CELL'] = self.cell
        os.environ['SGE_QMASTER_PORT'] = self.qmasterport
        os.environ['SGE_EXECD_PORT'] = self.execdport
        command = "/opt/sge6/bin/lx24-amd64/qconf -aattr hostgroup hostlist " + node.alias + " @allhosts >> /tmp/allhosts.out; "
        log.info("sge.addToAllhosts    command: %s", command)

        envars = self.exportEnvironmentVars()

        original_stdout = sys.stdout
        sys.stdout = NullDevice()  
        master.ssh.execute(envars + command)
        sys.stdout = original_stdout

    def addEnvarsToProfile(self, node):
        """
            Add environment variables (SGE_CELL, ports, etc.) to /etc/profile
        """
        envars = self.exportEnvironmentVars();
        log.debug("sge.addEnvarsToProfile    envars: echo '%s' >> /etc/profile", envars)
        node.ssh.execute("echo '" + envars + "' >> /etc/profile")

    def enableSchedulingInfo(self):
        """
            Enable job scheduling info output for 'qstat -j'
        """
        log.info("Enabling job scheduling info")

        envars = self.exportEnvironmentVars()
        log.debug(envars + "/opt/sge6/bin/lx24-amd64/qconf -ssconf")
        queue_template = subprocess.Popen(envars + "/opt/sge6/bin/lx24-amd64/qconf -ssconf", stdout=subprocess.PIPE, shell=True).stdout.read()
        log.debug("sge.CreateCell.enableSchedulingInfo    BEFORE queue_template: %s", queue_template)

        match = "schedd_job_info                   false"
        insert = "schedd_job_info                   true"
        queue_template = string.replace(queue_template, match, insert)
        log.debug("sge.CreateCell.enableSchedulingInfo    AFTER queue_template: %s", queue_template)

        pid = os.getpid()
        filename = "/tmp/queue-" + str(os.getpid()) + ".txt"
        queue_file = open(filename, 'w')
        print >> queue_file, queue_template
        queue_file.close()

        cmd = envars + "/opt/sge6/bin/lx24-amd64/qconf -Msconf " + filename
        log.debug(cmd)
        os.system(cmd)
        remove = "rm -fr " + filename
        log.debug(remove)
        os.system(remove)

    def addParallelEnvironment(self, master):
        """
            Add 'threaded' parallel environment
        """
        log.info("Adding 'threaded' parallel environment")

        sge_pe_template = """
        pe_name           threaded
        slots             %s
        user_lists        NONE
        xuser_lists       NONE
        start_proc_args   /bin/true
        stop_proc_args    /bin/true
        allocation_rule   $pe_slots
        control_slaves    TRUE
        job_is_first_task FALSE
        urgency_slots     min
        accounting_summary FALSE
        """
        pe_file = master.ssh.remote_file("/tmp/pe.txt")
        print >> pe_file, sge_pe_template % 99999
        pe_file.close()

        envars = self.exportEnvironmentVars()
        master.ssh.execute(envars + "/opt/sge6/bin/lx24-amd64/qconf -Ap %s &> /tmp/pe.out" % pe_file.name)
        master.ssh.execute(envars + '/opt/sge6/bin/lx24-amd64/qconf -mattr queue pe_list "threaded" all.q &> /tmp/pe2q.out')

    def addHeadToHosts(self, master):
        """
            Add head node to submit hosts and admin hosts lists
        """
        log.info("Adding head node to submit hosts and admin hosts lists")

        envars = self.exportEnvironmentVars()
        add_submit      = envars + '/opt/sge6/bin/lx24-amd64/qconf -as ' + self.head_ip
        add_admin       = envars + '/opt/sge6/bin/lx24-amd64/qconf -ah ' + self.head_ip

        log.debug("sge.addHeadToHosts    %s", add_submit)
        master.ssh.execute(add_submit)
        log.debug("sge.addHeadToHosts    %s", add_admin)
        master.ssh.execute(add_admin)

    def restartSge(self, node):
        """
            Restart SGE qmaster (master) and execd (master + nodes) daemons
        """
        log.info("Restarting SGE qmaster and execd daemons")

        envars = self.exportEnvironmentVars()
        killall = "/bin/ps aux | grep sgeadmin | cut -c9-14 | xargs -n1 -iPID /bin/kill -9 PID &> /dev/null"
        stop_execd      = envars + '/opt/sge6/bin/lx24-amd64/qconf -ke all'
        stop_qmaster    = envars + '/opt/sge6/bin/lx24-amd64/qconf -km'
        start_qmaster   = envars + '/opt/sge6/bin/lx24-amd64/sge_qmaster'
        start_execd     = envars + '/opt/sge6/bin/lx24-amd64/sge_execd'

        sleep = 1
        log.debug("sge.restartSge    Doing RESTART SGE: %s (%s)", node.alias, node.private_ip_address)

        log.info(killall)
        node.ssh.execute(killall, True, False, True)

        if node.alias == "master":
            time.sleep(float(sleep))
            log.debug("sge.restartSge    %s", start_qmaster)
            node.ssh.execute(start_qmaster)

        log.debug("sge.restartSge    %s", start_execd)
        node.ssh.execute(start_execd)

    def settingsCommand(self):
        target      =   self.root + "/" + self.cell + "/common"
        cmd     =   'cd ' + target + '; '
        cmd     +=  self.exportEnvironmentVars()
        cmd     +=  self.root + '/util/create_settings.sh ' + target
        log.debug("CreateCell.createSettings    cmd: %s", cmd)
        return cmd

    def createSettings(self, node):
        """
            Generate settings.sh file containing SGE_CELL, SGE_ROOT and port info
        """    
        log.info("Generating settings.sh file")
        log.debug("CreateCell.createSettings    CreateCell.createSettings(master)")
        cmd = self.settingsCommand()
        log.debug("CreateCell.createSettings    cmd: %s", cmd)
        node.ssh.execute(cmd)

    def exportEnvironmentVars(self):
        vars    =   'export SGE_ROOT=' + self.root + '; '
        vars    +=  'export SGE_CELL=' + self.cell + '; '
        vars    +=  'export SGE_QMASTER_PORT=' + self.qmasterport + '; '
        vars    +=  'export SGE_EXECD_PORT=' + self.execdport + '; '
        return vars

    def updateHeadIp(self):
        """
            Set hostname as head_ip (in case has changed due to reboot)
        """
        log.info("Updating hostname on head node")

        log.debug("CreateCell.updateHeadIp    self.head_long_ip: %s", self.head_long_ip)
        cmd = "hostname " + self.head_long_ip
        log.debug("CreateCell.updateHeadIp    cmd: %s", cmd)
        os.system(cmd)

    def updateMasterIp(self, master):
        """
            Replace 'master' with correct IP address in act_qmaster file
        """
        log.info("Updating act_qmaster file")
        target      =   self.root + "/" + self.cell
        act_qmaster =   target + "/common/act_qmaster"
        master_ip   =   master.private_ip_address

        log.debug("CreateCell.updateMasterIp    CreateCell.updateMasterIp(nodes)")
        log.debug("CreateCell.updateMasterIp    act_qmaster: %s", act_qmaster)
        log.debug("CreateCell.updateMasterIp    master_ip: %s", master_ip)
        cmd = "echo '" + master_ip + "' > " + act_qmaster
        log.debug("CreateCell.updateMasterIp    cmd: %s", cmd)
        os.system(cmd)

    def copyCellCommands(self):
        source      = self.root + "/default"
        target      = self.root + "/" + self.cell
        return (
            'mkdir ' + target + ' &> /dev/null', 
            'rsync -a ' + source + "/* " + target + " --exclude *tar.gz",
            'chown -R sgeadmin:sgeadmin ' + target
        )

    def copyCellOnHead(self):
        """
            Copy cell dir from default dir
        """
        log.info("Copying cell directory on head node")
        log.debug("CreateCell.copyCellOnHead    CreateCell.copyCellOnHead()")
        commands = self.copyCellCommands()
        log.debug("CreateCell.copyCellOnHead    commands: %s", commands)

        target      = self.root + "/" + self.cell
        log.debug("CreateCell.copyCell    target: %s", target)
        log.debug("CreateCell.copyCellOnHead    os.path.isdir(target): %s", os.path.isdir(target))
        if not os.path.isdir(target):
            for command in commands:
                log.info(command)
                os.system(command)

            ##### CREATE NEW settings.sh FILE
            os.system(self.settingsCommand)

    def copyCell(self, node):
        """
            Copy cell dir from default dir
        """
        log.info("Copying cell directory on %s", node.alias)
        log.debug("CreateCell.copyCell    CreateCell.copyCell(" + node.alias + ")")
        commands = self.copyCellCommands()
        log.debug("CreateCell.copyCell    commands: %s", commands)

        target      = self.root + "/" + self.cell
        log.debug("CreateCell.copyCell    target: %s", target)
        log.debug("CreateCell.copyCell    os.path.isdir(target): %s", os.path.isdir(target))
        #if not os.path.isdir(target):
        for command in commands:
            log.info(command)
            node.ssh.execute(command, True, False, True)

            ##### CREATE NEW settings.sh FILE
            node.ssh.execute(self.settingsCommand())

    def openSgePorts(self):
        """
            Open the particular SGE qmaster and execd daemon ports for this cluster
        """
        log.info("Opening SGE qmaster and execd ports")
        qmasterport = self.qmasterport
        execdport = self.execdport
        cluster = self.cell

        log.debug("CreateCell.openSgePorts    CreateCell.openSgePorts()")
        log.debug("CreateCell.openSgePorts    qmasterport; %s", qmasterport)
        log.debug("CreateCell.openSgePorts    execdport; %s", execdport)

        # HEAD NODE (I.E., NOT MASTER OR NODE)
        log.debug("CreateCell.openSgePorts    Opening ports")
        log.debug('ec2-authorize @sc-' + cluster + ' -p ' + execdport + ' -P tcp')
        os.system('ec2-authorize @sc-' + cluster + ' -p ' + execdport + ' -P tcp')
        log.debug('ec2-authorize @sc-' + cluster + ' -p ' + execdport + ' -P udp')
        os.system('ec2-authorize @sc-' + cluster + ' -p ' + execdport + ' -P udp')
        log.debug('ec2-authorize @sc-' + cluster + ' -p ' + qmasterport + ' -P tcp')
        os.system('ec2-authorize @sc-' + cluster + ' -p ' + qmasterport + ' -P tcp')
        log.debug('ec2-authorize @sc-' + cluster + ' -p ' + qmasterport + ' -P udp')
        os.system('ec2-authorize @sc-' + cluster + ' -p ' + qmasterport + ' -P udp')

    def on_add_node(self, node, nodes, master, user, user_shell, volumes):
        log.info("Doing 'on_add_node' for plugin: sge.CreateCell");
        log.info("Adding %s", node.alias)
        log.debug("CreateCell.on_add_node    CreateCell.on_add_node(self, node, nodes, master, user, user_shell, volumes)")
        log.debug("CreateCell.on_add_node    node.private_dns_name: %s" % node.private_dns_name)
        #### ADD ENVIRONMENT VARIABLES TO /etc/profile ON MASTER
        self.addEnvarsToProfile(node)

        ##### CREATE NEW CELL DIRECTORY ON HEAD AND MASTER
        self.copyCell(node);

        ##### RESTART SGE ON NODE
        self.restartSge(node)

        #### ADD NODE TO @allhosts GROUP
        self.addToAllhosts(node, master)

        log.info("Completed 'on_add_node' for plugin: sge.CreateCell");

    def on_remove_node(self, node, nodes, master, user, user_shell, volumes):
        log.info("Doing on_remove_node for plugin: sge.CreateCell")
        log.info("Removing %s " % node.alias)
        log.debug("CreateCell.on_remove_node    node.private_dns_name: %s" % node.private_dns_name)

