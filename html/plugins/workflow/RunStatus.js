dojo.provide("plugins.workflow.RunStatus");

// DISPLAY THE STATUS OF A WORKFLOW STAGE

// TITLE PANE
dojo.require("dijit.TitlePane");

// HAS A
dojo.require("plugins.core.ConfirmDialog");

// INHERITS
dojo.require("plugins.core.Common");

dojo.declare( "plugins.workflow.RunStatus",
	[ dijit._Widget, dijit._Templated, plugins.core.Common ],
{
//Path to the template of this widget. 
templatePath: dojo.moduleUrl("plugins", "workflow/templates/runstatus.html"),

// Calls dijit._Templated.widgetsInTemplate
widgetsInTemplate : true,

// OR USE @import IN HTML TEMPLATE
cssFiles : [ dojo.moduleUrl("plugins") + "/workflow/css/runstatus.css" ],

srcNodeRef: null,

timer: null,

// CORE WORKFLOW OBJECTS
core : null,

/////}
constructor : function(args) {

	// GET ARGS
	this.core = args.core;
	this.parentWidget = args.parentWidget;
	this.attachNode = args.attachNode;

	// LOAD CSS
	this.loadCSS();		
},

postMixInProperties: function() {
},

postCreate: function() {

	this.startup();
},

startup : function () {

	// SET UP THE ELEMENT OBJECTS AND THEIR VALUE FUNCTIONS
	this.inherited(arguments);


	// ADD TO TAB CONTAINER		
	this.attachNode.addChild(this.statusTab);
	//this.attachNode.selectChild(this.statusTab);	

	// START UP CONFIRM DIALOGUE
	this.setConfirmDialog();

	// INSTANTIATE NEW TIMER
	this.timer = new dojox.timing.Timer;
	this.timer.setInterval(15000);	// 1000ths OF A SECOND
},


getStatus : function () {
// KEEP POLLING SERVER FOR RUN STATUS UNTIL COMPLETE

	// SET MESSAGE
	this.notifier.innerHTML = "Getting run status...";

	if ( this.runner == null )
	{
	}

	var project		=	this.runner.project;
	var workflow	=	this.runner.workflow;
	var number		=	this.runner.number;
	var childNodes	=	this.runner.childNodes;

	// SANITY CHECKS
	if ( project == null )	return;
	if ( workflow == null )	return;
	if ( number == null )	return;
	if ( childNodes == null )	return;
	if ( ! childNodes || childNodes.length == 0 )
	{
		return;
	}

	// GET URL 
	var url = Agua.cgiUrl + "workflow.cgi";

	// GENERATE QUERY JSON FOR THIS WORKFLOW IN THIS PROJECT
	var query = new Object;
	query.username = Agua.cookie('username');
	query.sessionId = Agua.cookie('sessionId');
	query.project = project;
	query.workflow = workflow;
	query.number = number;
	query.mode = "getStatus";

	// RUN STATUS VARIABLES
	var completed = false;
	this.polling = true;

	// SET WORKFLOW OBJECT FOR CALL TO displayStatus
	var thisObject = this;

	this.timer.onTick = function()
	{

		var timer = this;

		// HANDLE RUN STATUS REPONSE FROM SERVER. POLLING 
		// WILL STOP WHEN THE LAST RUNNING ITEM HAS
		// COMPLETED, I.E., WHEN 'completed' == true

		if ( completed === true )
		{
			timer.stop();
		}

		// QUERY RUN STATUS ON SERVER
		dojo.xhrPut(
			{
				url: url,
				putData: dojo.toJson(query),
				handleAs: "json",
				//handleAs: "json-comment-optional",
				timeout : 5000,
				sync: false,
				handle: function (response)
				{
					if ( response == null )	return;


					// STORE THIS response
					this.statusResponse = response;

					// SET COMPLETED FLAG
					completed = true;

					// SET THE NODE CLASSES BASED ON STATUS

					// CHANGE CSS ON RUN NODES
					var startIndex = number - 1;
					for ( var i = startIndex; i < response.length; i++ )
					{
						var nodeClass = response[i].status;
						dojo.removeClass(childNodes[i], 'waiting');
						dojo.removeClass(childNodes[i], 'running');
						dojo.removeClass(childNodes[i], 'completed');
						dojo.addClass(childNodes[i], nodeClass);

						// UNSET COMPLETED FLAG IF ANY NODE IS NOT COMPLETED
						if ( nodeClass != "completed" )
						{
							completed = false;
						}
					}

					thisObject.showStatus(startIndex, response);
				}
			}
		);

	};	//	timer.onTick

	this.timer.start();
},


showStatus : function (startIndex, response) {
// POPULATE THE 'STATUS' PANE WITH RUN STATUS INFO

	console.clear();
	// SET TIME AT TOP OF STATUS PANEL
	this.setTime();

	// SHOW STAGES STATUS
	this.displayStagesStatus(response.stages);

	// SHOW CLUSTER STATUS
	this.displayClusterStatus(response.cluster);

	// SHOW QUEUE STATUS
	this.displayQueueStatus(response.queue);
},

setTime : function () {
	// SET DATE
	var timeNow = new Date();
	var timeString = timeNow.toString();
	//timeString.replace(/\s+GMT.+$/,'');
	timeString = timeString.match(/^(.+?)\s+GMT/)[1];
	this.notifier.innerHTML = timeString;			
},


displayClusterStatus : function (cluster) {

	// REMOVE ANY EXISTING STATUS TABLE
	while ( this.clusterStatusContainer.firstChild )
	{
		this.clusterStatusContainer.removeChild(this.clusterStatusContainer.firstChild);
	}

	this.clusterStatusContainer.innerHTML = "<PRE>" + cluster + "</PRE>";
},

displayQueueStatus : function (queue) {

	// REMOVE ANY EXISTING STATUS TABLE
	while ( this.queueStatusContainer.firstChild )
	{
		this.queueStatusContainer.removeChild(this.queueStatusContainer.firstChild);
	}

	this.queueStatusContainer.innerHTML = "<PRE>" + queue + "</PRE>";
},

displayStagesStatus : function (stages) {

	// REMOVE ANY EXISTING STATUS TABLE
	while ( this.stagesStatusContainer.firstChild )
	{
		this.stagesStatusContainer.removeChild(this.stagesStatusContainer.firstChild);
	}

	// BUILD TABLE
	var table = document.createElement('table');
	this.stagesStatusContainer.appendChild(table);

	// SET THE NODE CLASSES BASED ON STATUS
	for ( var i = startIndex; i < stages.length; i++ )
	{
		var tr = document.createElement('tr');
		table.appendChild(tr);
		dojo.addClass(table, 'displayStatusTable');

		var currentDate = new Date();
		var duration = '';
		// IF STARTED, GET THE DURATION
		if ( stages[i].started != null
			&& stages[i].started != "0000-00-00 00:00:00" )
		{
			// GET DURATION BY SUBSTRACTING started FROM completed OR now
			var startedDate = this.stringToDate(stages[i].started);
			var currentDate;
			if ( stages[i].completed == null
				|| stages[i].completed == ''
				|| stages[i].completed == "0000-00-00 00:00:00" )
			{
				currentDate = this.stringToDate(stages[i].now);
			}
			else
			{
				currentDate = this.stringToDate(stages[i].completed);
			}

			// CONVERT DIFFERENCE TO DURATION
			var seconds = currentDate - startedDate;
			duration = this.secondsToDuration(seconds);
		}
		stages[i].duration = duration;

		var statusRow = new plugins.workflow.RunStatusRow(stages[i]);

		var td = document.createElement('td');
		tr.appendChild(td);
		td.appendChild(statusRow.domNode);
	}
},

stringToDate : function (string) {
// CONVERT 2010-02-21 10:45:46 TO year, month, day, hour, minutes, seconds

	// CONVERT STRING TO ARRAY
	var array = string.match(/^(\d{4})\D+(\d+)\D+(\d+)\D+(\d+)\D+(\d+)\D+(\d+)/);

	// REMOVE FIRST MATCH (I.E., ALL OF MATCHED STRING)
	array.shift();

	// MONTH IS ZERO-INDEXED
	array[1]--;

	// GENERATE NEW DATE
	var date = new Date(array[0],array[1],array[2],array[3],array[4],array[5]);

	return date;
},


secondsToDuration : function (milliseconds) {
// CONVERT MILLISECONDS TO hours, mins AND secs 
	var duration = '';
	var remainder = 0;

	// IGNORE MILLISECONDS	
	remainder = milliseconds % 1000;
	//sComp = remainder.toString();
	//while (sComp.length<3)
	//	sComp = "0" + sComp;
	//duration = "." + sComp + " sec";
	milliseconds -= remainder;
	milliseconds /= 1000;

	// GET SECONDS	
	remainder = milliseconds % 60;
	if (remainder)
		duration = remainder.toString();
	else
		duration = "0";

	duration = duration + " sec";

	// Strip off last component
	milliseconds -= remainder;
	milliseconds /= 60;

	// GET HOURS	
	remainder = milliseconds % 60;
	duration = remainder.toString() + " min " + duration;

	// Strip off last component
	milliseconds -= remainder;
	milliseconds /= 60;

	// GET DAYS
	return milliseconds.toString() + " hours " + duration;
},


stopRun : function () {

	this.timer.stop();

	// SEND STOP SIGNAL TO SERVER
	this.sendStop();	
},	

pauseRun : function () {
// STOP AT THE CURRENT STAGE. TO RESTART FROM THIS STAGE, HIT 'START' BUTTON

	this.timer.stop();

	// SEND STOP SIGNAL TO SERVER
	this.sendStop();	

	// STORE THIS response
	var response = this.statusResponse;

	// SET THE NODE CLASSES BASED ON STATUS

	// CHANGE CSS ON RUN NODES
	var startIndex = number - 1;
	for ( var i = startIndex; i < response.length; i++ )
	{
		// SET this.runner.number TO FIRST RUNNING OR WAITING STAGE
		if ( response[i].status == "completed" )	continue;
		this.runner.number = (startIndex + 1);

		dojo.removeClass(childNodes[i], 'waiting');
		dojo.removeClass(childNodes[i], 'running');
		dojo.removeClass(childNodes[i], 'completed');
		dojo.addClass(childNodes[i], 'waiting');
		break;
	}
},

sendStop : function () {

	var project		=	this.runner.project;
	var workflow	=	this.runner.workflow;
	var number		=	this.runner.number;

	// SET TIMER CSS 
	dojo.removeClass(this.toggle, 'timerStopped');
	dojo.addClass(this.toggle, 'timerStarted');

	// SET MESSAGE
	this.notifier.innerHTML = "Running workflow...";

	// GET URL 
	var url = Agua.cgiUrl + "workflow.cgi";

	// GENERATE QUERY JSON FOR THIS WORKFLOW IN THIS PROJECT
	var query = new Object;
	query.username = Agua.cookie('username');
	query.sessionId = Agua.cookie('sessionId');
	query.project = project;
	query.workflow = workflow;
	query.mode = "stopWorkflow";
	query.number = number;

	dojo.xhrPut(
		{
			url: url,
			putData: dojo.toJson(query),
			handleAs: "json",
				//handleAs: "json-comment-optional",
			async: true,
			timeout: 100,
			handle: function(response){
				// DO NOTHING
			}
		}
	);

},


toggleTimer : function () {
// START TIME IF STOPPED OR STOP IT IF ITS RUNNING 

	if ( this.polling == true )	this.stopTimer();
	else this.startTimer();
},

stopTimer : function () {
// STOP POLLING THE SERVER FOR RUN STATUS

	if ( this.runner == null )	return;
	this.timer.stop();
	this.polling = false;

	dojo.removeClass(this.toggle, 'timerStarted');
	dojo.addClass(this.toggle, 'timerStopped');
},	

startTimer : function () {
// RESTART POLLING THE SERVER FOR RUN STATUS

	if ( this.runner == null )	return;
	this.getStatus();
	this.polling = true;

	dojo.removeClass(this.toggle, 'timerStopped');
	dojo.addClass(this.toggle, 'timerStarted');
},	


startRun : function (runner) {
// PROMPT FOR RESTART IF ALREADY RUNNING, OTHERWISE RUN WORKFLOW
	//console.dir(runner);

	// SELECT THIS TAB NODE
	this.attachNode.selectChild(this.statusTab);

	// REMOVE this.runner
	if ( this.runner != null )	delete this.runner;

	// SET this.runner
	this.runner = runner;

	var project		=	runner.project;
	var workflow	=	runner.workflow;
	var number		=	runner.number;


	// RUN THE WORKFLOW IF ITS NOT ALREADY RUNNING
	var isRunning = this.isRunning(project, workflow);

	// IF NOT STAGES ARE CURRENTLY RUNNING, RUN THE WORKFLOW
	if ( isRunning == false )	this.runWorkflow();

	// IF ITS ALREADY RUNNING, ASK FOR CONFIRMATION TO STOP AND RESTART
	if ( isRunning == true )
	{

		// CALLBACKS
		var noCallback = function (){
		};
		var Callback = function (){
		};
		var yesCallback = dojo.hitch(this, function()
			{
				this.stopRun();
				this.runWorkflow();
			}								
		);

		// SET TITLE AND MESSAGE
		var title = "Project '" + project + "' workflow '" + workflow + "' is already running";
		var message = "Do you want to stop and restart it?";

		// SHOW THE DIALOG
		this.loadConfirmDialog(title, message, yesCallback, noCallback);
	}

},

runWorkflow : function () {
// EXECUTE THE WORKFLOW

	var project		=	this.runner.project;
	var workflow	=	this.runner.workflow;
	var workflownumber	=	this.runner.workflownumber;
	var number		=	this.runner.number;
	var cluster		=	this.runner.cluster;

	// SET TIMER CSS 
	dojo.removeClass(this.toggle, 'timerStopped');
	dojo.addClass(this.toggle, 'timerStarted');

	// SET MESSAGE
	this.notifier.innerHTML = "Running workflow...";

	// GET URL 
	var url = Agua.cgiUrl + "workflow.cgi";

	// GENERATE QUERY JSON FOR THIS WORKFLOW IN THIS PROJECT
	var query = new Object;
	query.username = Agua.cookie('username');
	query.sessionId = Agua.cookie('sessionId');
	query.project = project;
	query.workflow = workflow;
	query.workflownumber = workflownumber;
	query.mode = "executeWorkflow";
	query.number = number;

	dojo.xhrPut(
		{
			url: url,
			putData: dojo.toJson(query),
			handleAs: "json",
				//handleAs: "json-comment-optional",
			async: true,
			timeout: 0,
			handle: function(response){
				// DO NOTHING
			}
		}
	);

	setTimeout(function(thisObj) { thisObj.getStatus(); }, 3000, this);
},

isRunning : function (project, workflow) {
// PERIODICALLY CHECK THE STATUS OF THE WORKFLOW

	// GET URL 
	var url = Agua.cgiUrl + "workflow.cgi";

	// GENERATE QUERY JSON FOR THIS WORKFLOW IN THIS PROJECT
	var query = new Object;
	query.username = Agua.cookie('username');
	query.sessionId = Agua.cookie('sessionId');
	query.project = project;
	query.workflow = workflow;
	query.mode = "getStatus";

	var isRunning = false;
	dojo.xhrPut(
		{
			url: url,
			putData: dojo.toJson(query),
			handleAs: "json",
				//handleAs: "json-comment-optional",
			sync: true,
			timeout: 100,
			handle: function(response){

				for ( var i = 0; i < response.length; i++ )
				{
					if ( response[i].status == "running" )
					{
						isRunning = true;
						break;
					}
				}
			}
		}
	);

	return isRunning;
},

setConfirmDialog : function () {
	var yesCallback = function (){};
	var noCallback = function (){};
	var title = "Dialog title";
	var message = "Dialog message";

	this.confirmDialog = new plugins.core.ConfirmDialog(
		{
			title 				:	title,
			message 			:	message,
			parentWidget 		:	this,
			yesCallback 		:	yesCallback,
			noCallback 			:	noCallback
		}			
	);
},

loadConfirmDialog : function (title, message, yesCallback, noCallback) {

	this.confirmDialog.load(
		{
			title 				:	title,
			message 			:	message,
			yesCallback 		:	yesCallback,
			noCallback 			:	noCallback
		}			
	);
}

});	// plugins.workflow.RunStatus


dojo.declare( "plugins.workflow.RunStatusRow",
[ dijit._Widget, dijit._Templated ],
{

//Path to the template of this widget. 
templatePath: dojo.moduleUrl("plugins", "workflow/templates/runstatusrow.html"),

// Calls dijit._Templated.widgetsInTemplate
widgetsInTemplate : true,

srcNodeRef: null,

//buildRendering: function(){
//	// we already have the srcNodeRef, so lets just
//	// keep it and use it as the domNode
//	this.domNode = this.srcNodeRef;
//
//	// call this._attachTemplateNodes to parse the template,
//	// which is actually just the srcnode
//	// this method is provided by the _Templated mixin
//	this._attachTemplateNodes(this.domNode);
//}

constructor : function(args) {	
},

postMixInProperties: function() {
},

postCreate : function() {

	this.startup();
},

startup : function () {

	this.inherited(arguments);

}



}); // plugins.workflow.RunStatusRow

