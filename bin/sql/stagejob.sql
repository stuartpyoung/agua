CREATE TABLE IF NOT EXISTS stagejob (

    username        VARCHAR(30) NOT NULL,
    project            VARCHAR(20) NOT NULL,
    workflow        VARCHAR(20) NOT NULL,

    name                VARCHAR(40) NOT NULL default '',    
    number              VARCHAR(10),

    type                VARCHAR(40),
    location            VARCHAR(255) NOT NULL default '',
    executor            VARCHAR(40) NOT NULL default '',
    cluster             VARCHAR(20)NOT NULL default '',    

    started             DATETIME,
    completed           DATETIME,
    workflowpid         INT(12),
    parentpid           INT(12),
    childpid            INT(12),
    status              VARCHAR(20),

    PRIMARY KEY  (username, project, workflow, name, number)
);
