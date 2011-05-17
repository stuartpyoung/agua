CREATE TABLE IF NOT EXISTS stage (

    username            VARCHAR(30) NOT NULL,
    owner               VARCHAR(30) NOT NULL,
    project             VARCHAR(20) NOT NULL,
    workflow            VARCHAR(20) NOT NULL,
    workflownumber      INT(12),

    name                VARCHAR(40) NOT NULL default '',    
    number              VARCHAR(10),

    type                VARCHAR(40),
    location            VARCHAR(255) NOT NULL default '',
    executor            VARCHAR(40) NOT NULL default '',
    cluster             VARCHAR(20)NOT NULL default '',
    submit              INT(1),

    stderrfile          varchar(255) default NULL,
    stdoutfile          varchar(255) default NULL,

    queued              DATETIME,
    started             DATETIME,
    completed           DATETIME,
    workflowpid         INT(12),
    stagepid            INT(12),
    stagejobid          INT(12),
    status              VARCHAR(20),

    stagename           VARCHAR(40),
    stagedescription    TEXT,
    stagenotes          TEXT,

    PRIMARY KEY  (username, project, workflow, number)
);