CREATE TABLE IF NOT EXISTS workflowparameter
(
    username        VARCHAR(30) NOT NULL,
    project         VARCHAR(20) NOT NULL,
    workflow        VARCHAR(20) NOT NULL,

    appname         VARCHAR(40),
    apptype         VARCHAR(40),

    name            VARCHAR(40) NOT NULL default '',
    category        VARCHAR(40) NOT NULL default '',
    type            VARCHAR(20) NOT NULL default '',
    argument        VARCHAR(40) NOT NULL default '',
    value           TEXT,
    discretion        VARCHAR(10) NOT NULL default '',
    format          VARCHAR(40),
    description     TEXT, 
    args            TEXT,
    params          TEXT,
    paramFunction   TEXT,

    PRIMARY KEY  (username, project, workflow, appname, name)
);
