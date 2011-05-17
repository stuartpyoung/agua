package Agua::Common::SGE;
use Moose::Role;
use Moose::Util::TypeConstraints;

=head2

	PACKAGE		Agua::Common::SGE

	PURPOSE

		CLUSTER METHODS FOR Agua::Common

=cut


use Data::Dumper;
use File::Path;

sub getEnvars {
	my $self		=	shift;
	my $username	=	shift;
	my $cluster		=	shift;


	my $query = qq{SELECT qmasterport
FROM clustervars
WHERE username = '$username'
AND cluster = '$cluster'};
	my $qmasterport = $self->dbobject()->query($query);
	print "Agua::Common::SGE::getEnvars    qmasterport not defined\n" and exit if not defined $qmasterport;

	my $execdport = $qmasterport + 1;
	my $sgeroot = $self->conf()->getKeyValue('cluster', 'SGEROOT');

	my $envars = {};
	$envars->{qmasterport}->{value} = $qmasterport;
	$envars->{execdport}->{value} 	= $execdport;
	$envars->{sgeroot}->{value} 	= $sgeroot;
	$envars->{sgecell}->{value} 	= $cluster;
	$envars->{qmasterport}->{export}= "export SGE_QMASTER_PORT=$qmasterport; ";
	$envars->{execdport}->{export} 	= "export SGE_EXECD_PORT=$execdport; ";
	$envars->{sgeroot}->{export} 	= "export SGE_ROOT=$sgeroot; ";
	$envars->{sgecell}->{export} 	= "export SGE_CELL=$cluster; ";
	$envars->{tostring} 	= "export SGE_QMASTER_PORT=$qmasterport; export SGE_EXECD_PORT=$execdport; export SGE_ROOT=$sgeroot; export SGE_CELL=$cluster";

	return $envars;
}

sub createQueuefile {
	my $self		=	shift;
	my ($queue, $queuefile, $slots) = @_;

	my ($queuefiledir) = $queuefile =~ /^(.+?)\/([^\/]+)$/;
	print "Agua::Common::SGE::createQueuefile    queuefiledir: $queuefiledir\n";
	File::Path::mkpath($queuefiledir) if not -d $queuefiledir;
	print "Agua::Common::SGE::createQueuefile    Can't create queuefiledir: $queuefiledir\n" and exit if not -d $queuefiledir;

	my $contents = $self->setQueuefileContents($queue, $slots);
	open(OUT, ">$queuefile") or die "Agua::Common::SGE::createQueuefile    Can't open queuefile: $queuefile\n";
	print OUT $contents;
	close(OUT) or die "Agua::Common::SGE::createQueuefile    Can't close queuefile: $queuefile\n";

	return $queuefile;
}

sub setQueuefileContents {
	my $self		=	shift;
	my $queue		=	shift;
	my $slots		=	shift;

	$slots = 1 if not defined $slots;
	print "Agua::Common::SGE::setQueuefileContents    queue: $queue\n";

	return qq{	
	qname                 $queue
	hostlist              \@allhosts
	seq_no                0
	load_thresholds       np_load_avg=20
	suspend_thresholds    NONE
	nsuspend              1
	suspend_interval      00:05:00
	priority              0
	min_cpu_interval      00:05:00
	processors            UNDEFINED
	qtype                 BATCH INTERACTIVE
	ckpt_list             NONE
	pe_list               make orte threaded
	rerun                 FALSE
	slots                 $slots
	tmpdir                /tmp
	shell                 /bin/sh
	prolog                NONE
	epilog                NONE
	shell_start_mode      posix_compliant
	starter_method        NONE
	suspend_method        NONE
	resume_method         NONE
	terminate_method      NONE
	notify                00:00:60
	owner_list            NONE
	user_lists            NONE
	xuser_lists           NONE
	subordinate_list      NONE
	complex_values        NONE
	projects              NONE
	xprojects             NONE
	calendar              NONE
	initial_state         default
	s_rt                  INFINITY
	h_rt                  INFINITY
	s_cpu                 INFINITY
	h_cpu                 INFINITY
	s_fsize               INFINITY
	h_fsize               INFINITY
	s_data                INFINITY
	h_data                INFINITY
	s_stack               INFINITY
	h_stack               INFINITY
	s_core                INFINITY
	h_core                INFINITY
	s_rss                 INFINITY
	h_rss                 INFINITY
	s_vmem                INFINITY
	h_vmem                INFINITY
};
}

sub getQueue {
	my $self		=	shift;
	my $queue		=	shift;
	print "Agua::Common::SGE::getQueue    queue not defined\n"
		and exit if not defined $queue;


	my $envars = $self->getEnvars($self->username(), $self->cluster());
	my $get_queue = "$envars->{tostring}; qconf -sq $queue";
	return `$get_queue`;
}

sub addQueue  {
	my $self		=	shift;
	my ($queue, $queuefile, $slots) = @_;

	$self->createQueuefile($queue, $queuefile);

	my $envars = $self->getEnvars($self->username(), $self->cluster());

	#### ADD THE QUEUE
	my $add = "$envars->{tostring}; qconf -Aq $queuefile";
	print "Agua::Common::SGE::addQueue    add: $add\n";
	print `$add`;
}

sub modifyQueue  {
	my $self		=	shift;
	my ($queue, $queuefile) = @_;

	print "Agua::Common::SGE::removeQueue    Can't find queuefile: $queuefile\n" and exit if not -f $queuefile;

	#### ADD THE QUEUE
	my $modify = "qconf -Mq $queuefile";
	print "Agua::Common::SGE::modifyQueue    modify: $modify\n";
	print `$modify`;
}


sub removeQueue {
	my $self		=	shift;
	my ($queue, $queuefile) = @_;
print "Agua::Common::SGE::removeQueue    Agua::Common::SGE::removeQueue(queue, queuefile)\n";
print "Agua::Common::SGE::removeQueue    queue: $queue\n";
print "Agua::Common::SGE::removeQueue    queuefile: $queuefile\n";

	#### DELETE THE QUEUE USING qconf
	my $delete = "qconf -dq $queue";
	print "Agua::Common::SGE::removeQueue    delete: $delete\n";
	print `$delete`;	

	#### REMOVE THE QUEUEFILE
	my $remove = "rm -fr $queuefile";
	print "Agua::Common::SGE::removeQueue    remove: $remove\n";
	print `$remove`;	
}


sub getHostlist {
	my $self	=	shift;

	print "Agua::Common::SGE::getHostlist    Agua::Common::SGE::getHostlist()\n";

	my $command = "qconf -sel";
	my $hosts = `$command`;
	print "Agua::Common::SGE::getHostlist    hosts: @$hosts\n";

	my $hostlist = join ",", @$hosts;
	print "Agua::Common::SGE::getHostlist    hostlist: $hostlist\n";

	return $hostlist;
}

sub addPE {
	my $self		=	shift;
	my ($pe, $pefile, $slots) = @_;
	print "Agua::Common::SGE::setPE    pefile: $pefile\n"; 

	my ($pefiledir) = $pefile =~ /^(.+?)\/([^\/]+)$/;
	print "Agua::Common::SGE::removeQueue    pefiledir: $pefiledir\n";
	File::Path::mkpath($pefiledir) if not -d $pefiledir;
	print "Agua::Common::SGE::removeQueue    Can't create pefiledir: $pefiledir\n" and exit if not -d $pefiledir;

	#### BACKUP IF EXISTS
	$self->backupFile($pefile) if -f $pefile;

	open(OUT, ">$pefile") or die "Agua::Common::SGE::addPE    Can't open pefile: $pefile\n";
	my $contents = $self->setPEFileContents($pe, $slots);
	print "Agua::Common::SGE::addPE    contents: $contents\n";
	print OUT $contents;
	close(OUT) or die "Agua::Common::SGE::addPE    Can't close pefile: $pefile\n";
	print "Agua::Common::SGE::addPE    pefile: $pefile\n";

	#### ADD THE QUEUE
	my $add = "qconf -Ap $pefile";
	print "Agua::Common::SGE::addPE    add: $add\n";
	print `$add`;	

} ####	addPE

sub setPEFileContents {
	my $self		=	shift;
	my $pe			=	shift;
	my $slots		=	shift;

	# OPTIONAL ENTRIES:
    # slots              AUTO

	return qq{	
    pe_name            $pe
    user_lists         arusers
    slots              $slots
    xuser_lists        NONE
    start_proc_args    /bin/true
    stop_proc_args     /bin/true
    allocation_rule    \$pe_slots
    control_slaves     TRUE
    job_is_first_task  FALSE
    urgency_slots      min
    accounting_summary TRUE\n};
}

sub addPEToQueue {
	my $self		=	shift;
	my ($pe, $queue, $queuefile) = @_;
	print "Agua::Common::SGE::addPEToQueue    pe: $pe\n";
	print "Agua::Common::SGE::addPEToQueue    queue: $queue\n";
	print "Agua::Common::SGE::addPEToQueue    queuefile: $queuefile\n";

	#### ADD threaded PE TO default QUEUE'S pe_list:
	my $contents = $self->addPEQueuefileContents($pe, $queue, $queuefile);

	#### PRINT MODIFIED CONTENT TO QUEUE CONF FILE
	open(OUT, ">$queuefile") or die "Monitor::SGE::addPEToQueue    Can't write to queuefile: $queuefile\n";
	print OUT $contents;
	close(OUT) or die "Monitor::SGE::addPEToQueue    Can't close queuefile: $queuefile\n";

	$self->modifyQueue($queue, $queuefile);
}

sub addPEQueuefileContents {
	my $self		=	shift;
	my ($pe, $queue, $queuefile) = @_;
	print "Agua::Common::SGE::addPEQueuefileContents    pe: $pe\n";
	print "Agua::Common::SGE::addPEQueuefileContents    queuefile: $queuefile\n";

	my $contents;
	$contents = $self->fileContents($queuefile) if -f $queuefile and not -z $queuefile;
	$contents = $self->setQueuefileContents($queue)
		if not -f $queuefile or -z $queuefile or $contents =~ /^\s*$/;

	#### BACKUP QUEUEFILE IF ITS PRESENT AND NON-EMPTY
	$self->backupFile($queuefile) if -f $queuefile and not -z $queuefile and $contents !~ /^\s*$/;

	print "Agua::Common::SGE::addPEQueuefileContents    contents: $contents\n";

	if ( $contents =~ /^(.+)\n(\s*pe_list\s+)([^\n]+)\n(.+)$/ms )
	{
		$contents = $1 . "\n". $2 . $3 . " " . $pe . "\n" . $4 ;
	}

	print "Agua::Common::SGE::addPEQueuefileContents    contents: $contents\n";

	return $contents;	
}

sub fileContents {
	my $self		=	shift;
	my $file		=	shift;

	open(FILE, $file) or die "Agua::Common::SGE::addPEQueuefileContents    Can't open file: $file\n";
	my $temp = $/;
	$/ = undef;
	my $contents = 	<FILE>;
	close(FILE);
	$/ = $temp;
}
sub removePE {
	my $self		=	shift;
	my ($pe, $pefile) = @_;

	print "Agua::Common::SGE::removePE    Agua::Common::SGE::removePE(pe, pefile)\n";
	print "Agua::Common::SGE::removePE    pe: $pe\n";
	print "Agua::Common::SGE::removePE    pefile: $pefile\n";

	#### DELETE THE QUEUE USING qconf
	my $delete = "qconf -dq $pe";
	print "Agua::Common::SGE::removePE    delete: $delete\n";
	print `$delete`;	

	#### REMOVE THE QUEUEFILE
	my $remove = "rm -fr $pefile";
	print "Agua::Common::SGE::removePE    remove: $remove\n";
	print `$remove`;	
}


sub backupFile {
	my $self		=	shift;
	my $filename 	=	shift;

	#### BACKUP FILE
	my $counter = 1;
	my $backupfile = "$filename.$counter";
	while ( -f $backupfile )
	{
		$counter++;
		$backupfile = "$filename.$counter";
	}
	print "AWS::addToFstab    backupfile: $backupfile\n";
	`cp $filename $backupfile`;

	print "AWS::addToFstab    backupfile not created\n" and exit if not -f $backupfile;
}



sub getQueuefile {
	my $self		=	shift;
	my $queue		=	shift;
	print "Agua::Common::Cluster::getQueuefile    queue: $queue\n";

	my $adminuser = $self->conf()->getKeyValue("agua", 'ADMINUSER');
	$adminuser = $self->conf()->getKeyValue('agua', "ADMINUSER") if not defined $adminuser;
	print "Agua::Common::SGE::setEnv    sgeroot not defined\n" and exit if not defined $adminuser;
	print "Agua::Common::Cluster::getQueuefile    admin: $adminuser\n";
	my $fileroot = $self->getFileroot($adminuser);
	my $queuefile = "$fileroot/.sge/conf/$queue.conf";
	print "Agua::Common::Cluster::getQueuefile    queuefile: $queuefile\n";

	return $queuefile;
}


sub sgebinCommand {
	my $self		=	shift;

	my $username 	=	$self->username();
	my $cluster 	=	$self->cluster();

	#### SET COMMAND WITH ENVIRONMENT VARIABLES
 	my $envars = $self->getEnvars($username, $cluster);
	my $sgebin = $self->conf()->getKeyValue("cluster", 'SGEBIN');

    my $command;	
	$command .= "$envars->{tostring};" if defined $envars and defined $envars->{tostring};
	$command .= "$sgebin";

	return $command;
}

sub queueStatus {
	my $self		=	shift;

	my $username 	=	$self->username();
	my $cluster 	=	$self->cluster();

    my $command = $self->sgebinCommand($username, $cluster) . "/qstat -f";

    my $output = `$command`;

#	my $output = qq{queuename                      qtype resv/used/tot. load_avg arch          states
#---------------------------------------------------------------------------------
#Project1-Workflow1\@master      BIP   0/1/1          1.13     lx24-amd64    
#     37 0.55500 test       root         t     05/17/2011 03:43:54     1        
#---------------------------------------------------------------------------------
#Project1-Workflow1\@node001     BIP   0/1/1          1.29     lx24-amd64    
#     38 0.55500 test       root         t     05/17/2011 03:43:54     1        
#---------------------------------------------------------------------------------
#all.q\@master                   BIP   0/1/1          1.13     lx24-amd64    
#     35 0.55500 test       root         r     05/17/2011 03:43:54     1        
#---------------------------------------------------------------------------------
#all.q\@node001                  BIP   0/1/1          1.29     lx24-amd64    
#     36 0.55500 test       root         r     05/17/2011 03:43:54     1        
#
#############################################################################
# - PENDING JOBS - PENDING JOBS - PENDING JOBS - PENDING JOBS - PENDING JOBS
#############################################################################
#     39 0.55500 test       root         qw    05/17/2011 03:42:32     1        
#     40 0.55500 test       root         qw    05/17/2011 03:42:32     1        
#     41 0.55500 test       root         qw    05/17/2011 03:42:32     1        
#     42 0.55500 test       root         qw    05/17/2011 03:42:33     1        
#     43 0.55500 test       root         qw    05/17/2011 03:42:33     1        
#     44 0.55500 test       root         qw    05/17/2011 03:42:33     1        
#     45 0.55500 test       root         qw    05/17/2011 03:42:33     1        
#     46 0.55500 test       root         qw    05/17/2011 03:42:33     1        
#     47 0.55500 test       root         qw    05/17/2011 03:42:33     1        
#     48 0.55500 test       root         qw    05/17/2011 03:42:33     1        
#     49 0.55500 test       root         qw    05/17/2011 03:42:33     1        
#     50 0.55500 test       root         qw    05/17/2011 03:42:33     1        
#     51 0.55500 test       root         qw    05/17/2011 03:42:33     1        
#     52 0.55500 test       root         qw    05/17/2011 03:42:33     1        
#     53 0.55500 test       root         qw    05/17/2011 03:42:33     1        
#     54 0.55500 test       root         qw    05/17/2011 03:42:33     1        
#     55 0.55500 test       root         qw    05/17/2011 03:43:46     1        
#     56 0.55500 test       root         qw    05/17/2011 03:43:46     1        
#     57 0.55500 test       root         qw    05/17/2011 03:43:46     1        
#     58 0.55500 test       root         qw    05/17/2011 03:43:46     1        
#     59 0.55500 test       root         qw    05/17/2011 03:43:46     1        
#     60 0.55500 test       root         qw    05/17/2011 03:43:46     1        
#     61 0.55500 test       root         qw    05/17/2011 03:43:46     1        
#     62 0.55500 test       root         qw    05/17/2011 03:43:46     1        
#     63 0.55500 test       root         qw    05/17/2011 03:43:46     1        
#     64 0.55500 test       root         qw    05/17/2011 03:43:46     1        
#     65 0.55500 test       root         qw    05/17/2011 03:43:46     1        
#     66 0.55500 test       root         qw    05/17/2011 03:43:46     1        
#     67 0.55500 test       root         qw    05/17/2011 03:43:46     1        
#     68 0.55500 test       root         qw    05/17/2011 03:43:46     1        
#     69 0.55500 test       root         qw    05/17/2011 03:43:46     1        
#     70 0.55500 test       root         qw    05/17/2011 03:43:46     1        
#     71 0.55500 test       root         qw    05/17/2011 03:43:46     1        
#     72 0.55500 test       root         qw    05/17/2011 03:43:46     1        
#     73 0.55500 test       root         qw    05/17/2011 03:43:46     1        
#     74 0.55500 test       root         qw    05/17/2011 03:43:46     1        
#     75 0.55500 test       root         qw    05/17/2011 03:43:46     1        
#     76 0.55500 test       root         qw    05/17/2011 03:43:46     1        
#     77 0.55500 test       root         qw    05/17/2011 03:43:46     1        
#     78 0.55500 test       root         qw    05/17/2011 03:43:46     1        
#     79 0.55500 test       root         qw    05/17/2011 03:43:46     1        
#     80 0.55500 test       root         qw    05/17/2011 03:43:46     1        
#     81 0.55500 test       root         qw    05/17/2011 03:43:46     1        
#     82 0.55500 test       root         qw    05/17/2011 03:43:46     1        
#     83 0.55500 test       root         qw    05/17/2011 03:43:46     1        
#     84 0.55500 test       root         qw    05/17/2011 03:43:46     1        
#     85 0.55500 test       root         qw    05/17/2011 03:43:46     1        
#     86 0.55500 test       root         qw    05/17/2011 03:43:46     1        
#     87 0.55500 test       root         qw    05/17/2011 03:43:47     1        
#     88 0.55500 test       root         qw    05/17/2011 03:43:47     1        
#     89 0.55500 test       root         qw    05/17/2011 03:43:47     1        
#     90 0.55500 test       root         qw    05/17/2011 03:43:47     1        
#     91 0.55500 test       root         qw    05/17/2011 03:43:47     1        
#     92 0.55500 test       root         qw    05/17/2011 03:43:47     1        
#     93 0.55500 test       root         qw    05/17/2011 03:43:47     1        
#     94 0.55500 test       root         qw    05/17/2011 03:43:47     1   
#};

    my ($nodes, $jobs) = $output =~ /^(.+)\n#{10,}\n.+?PENDING JOBS.+?\n#{10,}\n(.+)$/ms;

    my @lines = split "\n", $jobs;
    my $jobcounts = {};
    foreach my $line ( @lines )
    {
        my ($jobname) = $line =~ /^\s*\S+\s+\S+\s+(\S+)/;
        if ( exists $jobcounts->{$jobname} )
        {
            $jobcounts->{$jobname}++;
        }
        else
        {
            $jobcounts->{$jobname} = 1;
        }
    }

    my $summary = qq{=================================================================================
$nodes
PENDING JOBS:\n};

    foreach my $key ( sort keys %$jobcounts )
    {
        $summary .= "$key\t$jobcounts->{$key} jobs\n";
    }

	return $summary;
}

1;