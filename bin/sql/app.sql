CREATE TABLE IF NOT EXISTS app
(
    owner           VARCHAR(30) NOT NULL,    
    name            VARCHAR(40) NOT NULL,
    type            VARCHAR(40) NOT NULL,    
    location        VARCHAR(255) NOT NULL default '',
    localonly       INT(1) NOT NULL default 0,
    executor        VARCHAR(40) NOT NULL default '',
    version         VARCHAR(20)NOT NULL default '',

    description     TEXT,
    notes           TEXT,

    PRIMARY KEY  (owner, name)
);
