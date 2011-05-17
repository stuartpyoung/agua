CREATE TABLE IF NOT EXISTS monitor
(
	processid	VARCHAR(20) NOT NULL,
	command     TEXT,
	outputdir   TEXT,
	datetime	DATETIME NOT NULL,
	status		TEXT,
	PRIMARY KEY (processid, datetime)
);
