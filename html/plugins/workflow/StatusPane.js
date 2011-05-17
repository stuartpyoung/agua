dojo.provide("plugins.workflow.RunStatus");

// DISPLAY THE STATUS OF A WORKFLOW STAGE

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

// CORE WORKFLOW OBJECTS
core : null,

/////}
// constructor
// Any initialization code would go here in the constructor.
// plugins.report.Template and its superclasses dijit._Widget and
// dijit._Templated do not have parameters in their constructors, so
// there wouldn't be any multiple-inheritance complications
// if you were to include some paramters here.
constructor : function(args)
{		

	this.core = args.core;

	// GET ARGS
	this.parentWidget = args.parentWidget;
	this.attachNode = args.attachNode;
},

// postMixInProperties
//Inherited from dijit._Widget and called just before template
//instantiation in buildRendering. This method is especially useful
//for manipulating the template before it becomes visible.
postMixInProperties: function()
{
},

//postCreate
//You can override this method to manipulate widget once it is
//placed in the UI, but be warned that any child widgets contained
//in it may not be ready yet.        
postCreate: function()
{
},

// startup
//
startup : function ()
{

	// SET UP THE ELEMENT OBJECTS AND THEIR VALUE FUNCTIONS
	this.inherited(arguments);


	// ADD TO TAB CONTAINER		
	this.attachNode.addChild(this.mainTab);

},


// KEEP POLLING SERVER FOR RUN STATUS UNTIL COMPLETE
getStatus : function ()
{

	// RUN STARTED FROM CONTEXT MENU ON STAGE NODE
	if ( this.runner != null )
	{
		var project		=	this.runner.project;
		var workflow	=	this.runner.workflow;
		var number		=	this.runner.number;
		var childNodes	=	this.runner.childNodes;
	}
	// OR, RUN STARTED BY CLICK ON 'RUN' BUTTON
	else
	{
		var project		=	this.project;
		var workflow	=	this.workflow;
		var number		=	1;
		var childNodes	=	this.dropTarget.getAllNodes();
	}


	// RETURN IF NO CHILDNODES, I.E., NO APPLICATIONS IN THE WORKFLOW
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


	// DEBUG ONLY:
	////// MAYBE EXTEND INTERVAL INSTEAD OF LIMITING NUMBER OF TRIES
	////var limitTries = 1000;

	// RUN STATUS VARIABLES
	var counter = 0;
	var completed = false;
	this.polling = true;

	// INSTANTIATE NEW TIMER
	this.timer = new dojox.timing.Timer;
	this.timer.setInterval(60000);	// 1000ths OF A SECOND

	// SET WORKFLOW OBJECT FOR CALL TO displayStatus
	var thisObject = this;

	this.timer.onTick = function()
	{
		counter++;

		var timer = this;

		//if ( counter >= limitTries )
		//{
		//	timer.stop();
		//}


		// HANDLE RUN STATUS REPONSE FROM SERVER.
		// THE POLLING WILL STOP WHEN THE LAST
		// RUNNING ITEM HAS COMPLETED, I.E., WHEN
		// 'completed' == true
		var handleStatus = function(response)
		{
			if  ( ! response )
			{
				return;
			}


			// SET COMPLETED FLAG
			completed = true;

			// SET THE NODE CLASSES BASED ON STATUS
			for ( var i = start; i < response.length; i++ )
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

			thisObject.showStatus(start, response);
		};	//	handleStatus


		timer.stop();


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
				sync: false,
				handle: handleStatus
			}
		);

	};	//	timer.onTick

	this.timer.start();
},


// POPULATE THE 'STATUS' PANE WITH RUN STATUS INFO
showStatus : function (start, response)
{


	response = [{'childpid':'234835','number':'1','status':'','project':'Project1','workflowpid':null,'completed':null,'owner':'admin','cluster':'','location':'apps/clusterELAND.pl','executor':'perl','name':'ELAND','username':'admin','workflow':'Workflow1','parentpid':'5889','type':'alignment','started':'2010-02-18 12:13:03'},{'childpid':'234836','number':'2','status':'0','project':'Project1','workflowpid':null,'completed':null,'owner':'admin','cluster':'PBS','location':'apps/MAQ.pl','executor':'/usr/bin/perl','name':'MAQ','username':'admin','workflow':'Workflow1','parentpid':'5889','type':'pipeline','started':null},{'childpid':'234835','number':'3','status':'','project':'Project1','workflowpid':null,'completed':null,'owner':'admin','cluster':'','location':'apps/filterSNP.pl','executor':'perl','name':'filterSNP','username':'admin','workflow':'Workflow1','parentpid':'5889','type':'report','started':null}];


	// REMOVE ANY EXISTING STATUS TABLE
	while ( this.statusContainer.firstChild )
	{
		this.statusContainer.removeChild(this.statusContainer.firstChild);
	}

	// BUILD TABLE
	var table = document.createElement('table');
	this.statusContainer.appendChild(table);

	// SET THE NODE CLASSES BASED ON STATUS
	for ( var responseCounter = start; responseCounter < response.length; responseCounter++ )
	{
		var tr = document.createElement('tr');
		table.appendChild(tr);
		dojo.addClass(table, 'displayStatusTable');

		var statusRow = new plugins.workflow.RunStatusRow(
			{
				number: response[responseCounter].number,
				name: response[responseCounter].name,
				status: response[responseCounter].status,
				started: response[responseCounter].started,
				completed: response[responseCounter].completed
			}
		);

		var td = document.createElement('td');
		tr.appendChild(td);
		td.appendChild(statusRow.domNode);
	}
},


stopRun : function ()
{

},	


// EXECUTE THE RUN
startRun : function (runner)
{

	var project		=	runner.project;
	var workflow	=	runner.workflow;
	var number		=	runner.number;



	var running = this.checkRunning(project, workflow);


	if ( running.status != null && running.status == "running" )
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
	query.mode = "executeWorkflow";
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

	setTimeout(function(thisObj) { thisObj.getStatus(); }, 5000, this);
},



// PERIODICALLY CHECK THE STATUS OF THE WORKFLOW
checkRunning : function (project, workflow)
{

	// GET URL 
	var url = Agua.cgiUrl + "workflow.cgi";

	// GENERATE QUERY JSON FOR THIS WORKFLOW IN THIS PROJECT
	var query = new Object;
	query.username = Agua.cookie('username');
	query.sessionId = Agua.cookie('sessionId');
	query.project = project;
	query.workflow = workflow;
	query.mode = "getStatus";

	var running;
	dojo.xhrPut(
		{
			url: url,
			putData: dojo.toJson(query),
			handleAs: "json",
				//handleAs: "json-comment-optional",
			sync: true,
			timeout: 100,
			handle: function(response){

				running = response;

			}
		}
	);

	return running;
}

});	// plugins.workflow.RunStatus





dojo.declare( "plugins.workflow.RunStatusRow",
[ dijit._Widget, dijit._Templated ],
{

//Path to the template of this widget. 
templatePath: dojo.moduleUrl("plugins", "workflow/templates/historypanerow.html"),

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

constructor : function(args)
{		
},

// postMixInProperties
//Inherited from dijit._Widget and called just before template
//instantiation in buildRendering. This method is especially useful
//for manipulating the template before it becomes visible.
postMixInProperties: function()
{
},

postCreate : function()
{

	this.startup();
},

startup : function ()
{

	this.inherited(arguments);


}



});
