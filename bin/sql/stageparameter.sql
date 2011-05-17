CREATE TABLE IF NOT EXISTS stageparameter
(
    username        VARCHAR(30) NOT NULL,
    owner           VARCHAR(30) NOT NULL,
    project         VARCHAR(20) NOT NULL,
    workflow        VARCHAR(20) NOT NULL,
    workflownumber  INT(12),

    appname         VARCHAR(40),
    appnumber       VARCHAR(10),
    name            VARCHAR(40) NOT NULL default '',

    paramtype       VARCHAR(40) NOT NULL default '',
    category        VARCHAR(40) NOT NULL default '',
    valuetype       VARCHAR(20) NOT NULL default '',
    argument        VARCHAR(40) NOT NULL default '',
    value           TEXT,
    discretion      VARCHAR(10) NOT NULL default '',
    format          VARCHAR(40),
    description     TEXT, 
    args            TEXT,
    params          TEXT,
    paramFunction   TEXT,

    chained         INT(1),

    PRIMARY KEY  (username, project, workflow, appname, appnumber, name)
);
