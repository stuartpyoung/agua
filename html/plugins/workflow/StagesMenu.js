
dojo.provide("plugins.workflow.StagesMenu");

// PROVIDE A POPUP CONTEXT MENU AND IMPLEMENT ITS 
// ONCLICK FUNCTIONS

//Open a context menu	On Windows: shift-f10 or the Windows context menu key On Firefox on the Macintosh: ctrl-space. On Safari 4 on Mac: VO+shift+m (VO is usually control+opton)
//Navigate menu items	Up and down arrow keys
//Activate a menu item	Spacebar or enter
//Open a submenu	Spacebar, enter, or right arrow
//Close a context menu or submenu	Esc or left arrow
//Close a context menu and all open submenus	Tab

//dojo.require("dijit.dijit"); // optimize: load dijit layer
dojo.require("dojo.parser");

// HAS A
dojo.require("plugins.menu.Menu");

// INHERITS
dojo.require("plugins.core.Common");


dojo.declare(
    "plugins.workflow.StagesMenu",
	[ dijit._Widget, dijit._Templated, plugins.core.Common ],
{
//Path to the template of this widget. 
templatePath: dojo.moduleUrl("plugins", "workflow/templates/stagesmenu.html"),

// Calls dijit._Templated.widgetsInTemplate
widgetsInTemplate : true,

//addingApp STATE2
addingApp : false,

// OR USE @import IN HTML TEMPLATE
cssFiles : [ dojo.moduleUrl("plugins") + "/workflow/css/stagesmenu.css" ],

// PARENT WIDGET
parentWidget : null,

// CORE WORKFLOW OBJECTS
core : null,

////////}

constructor : function(args) {
	this.core = args.core;

	// GET INFO FROM ARGS
	this.parentWidget = args.parentWidget;

	// LOAD SORIA AND FILEPICKER CSS
	this.loadCSS();		
},

postCreate : function() {

	this.startup();
},

////////}

startup : function () {

	// COMPLETE CONSTRUCTION OF OBJECT
	this.inherited(arguments);	 

	// SET DRAG APP - LIST OF APPS
	this.setMenu();
},

// CONNECT LISTENERS FOR MENU
setMenu : function () {

	dojo.connect(this.removeNode, "onClick", dojo.hitch(this, function(event)
	{
		this.remove(event);
		event.stopPropagation();
	}));

	dojo.connect(this.runNode, "onClick", dojo.hitch(this, function(event)
	{
		this.run(event);
	}));
},

// BIND THE MENU TO A NODE
bind : function (node) {

	if ( node == null )
	{
	}
	return this.menu.bindDomNode(node);	
},

// CHAIN THE INPUTS AND OUTPUTS OF THIS STAGE TO THE PARAMETER VALUES
// OF THE PRECEDING STAGE
chain : function () {

	// REM: WE ARE NOT INTERESTED IN event.target 
	// BECAUSE ITS THE CLICKED MENU NODE. WE WANT
	// THE NODE UNDERNEATH
	var node = this.menu.currentTarget;
	var widget = dijit.getEnclosingWidget(node);
	var application = widget.getApplication();

	// PREPARE STAGE OBJECT
	var stageObject = {
		project: this.parentWidget.projectCombo.getValue(),
		workflow: this.parentWidget.workflowCombo.getValue(),
		owner: application.owner,
		appname: application.name,
		appnumber: application.number,
		name: application.name,
		number: application.number
	};


	// CHANGE THE STAGE PARAMETERS FOR THIS APPLICATION
	// IF THE args FIELD IS NOT NULL (ALSO params AND paramFunction)
	var force = true;
	this.parentWidget.workflowIO.chainStage(stageObject, force);

	// SET INFO PANE FOR DROPPED NODE

	var actualNode = widget.domNode.parentNode;
	//node.parentWidget = widget;
	this.core.stages.loadParametersPane(actualNode);
},

// REMOVE THE STAGE FROM THE WORKFLOW
remove : function (event) {

	// REM: WE ARE NOT INTERESTED IN event.target 
	// BECAUSE ITS THE CLICKED MENU NODE. WE WANT
	// THE NODE UNDERNEATH
	var node = this.menu.currentTarget;
	var widget = dijit.getEnclosingWidget(node);
	var application = widget.getApplication();

	if ( widget == null )	return;

	// PHYSICALLY REMOVE THE CLICKED NODE
	var itemNode = widget.domNode.parentNode;
	this.parentWidget.dropTarget.delItem(itemNode.id);
	dojo.destroy(itemNode);

	// SET username
	var username = Agua.cookie(username);

	// REMOVE THE CLICKED STAGE FROM THE WORKFLOW
	var stageObject = {
		username: username,
		project: this.parentWidget.projectCombo.getValue(),
		workflow: this.parentWidget.workflowCombo.getValue(),
		name: application.name,
		owner: application.owner,
		number: application.number,
		type: application.type
	};

	// REMOVE STAGE IN AGUA ON CLIENT AND ON REMOTE SERVER
	setTimeout(function(thisObj) {
			Agua.spliceStage(stageObject);
		},
		100,
		this
	);


	// UPDATE ANY NODES COMING AFTER THE INSERTION POINT OF THE NEW NODE
	// NB: THE SERVER SIDE UPDATES ARE DONE AUTOMATICALLY
	var childNodes = this.parentWidget.dropTarget.getAllNodes();
	for ( var i = application.number; i < childNodes.length; i++ )
	{
		childNodes[i].application.number = (i + 1).toString();
	}

	// RESETTING number IN ALL CHILDNODES
	for ( var i = 0; i < childNodes.length; i++ )
	{
		var node = childNodes[i];
		////console.dir(node);

		// GET WIDGET
		var widget = dijit.byNode(node.firstChild);
		//var widget = node.parentWidget;
		if ( widget == null )
		{
			widget = dijit.getEnclosingWidget(childNodes[i]);
		}
		node.application.number = (i + 1).toString();
		node.application.appnumber = (i + 1).toString();

		widget.setNumber(node.application.number);

	}


	// DO INFOPANE
	if ( childNodes.length )
	{
		this.parentWidget.loadParametersPane(childNodes[0]);
	}
	else
	{
		this.parentWidget.clearParameters();
	}

	// CALL Stages TO CHECK VALID STAGES		
	this.parentWidget.updateValidity();

},	//	remove


// ADD PROGRAMMATIC CONTEXT MENU
run : function () {

	var node = this.menu.currentTarget;
	var widget = dijit.getEnclosingWidget(node);
	var application = widget.getApplication();
	//var application = node.application;

	if ( application == null )
	{
		var widget = node.parentWidget;
		application = widget.getApplication();
	}

	var childNodes = this.parentWidget.dropTarget.getAllNodes();

	// INSTANTIATE RUNNER CLASS
	var runner = new Object;
	runner.project = this.parentWidget.getProject();
	runner.workflow = this.parentWidget.getWorkflow();
	runner.cluster = this.parentWidget.getCluster();
	runner.workflownumber = this.parentWidget.getWorkflowNumber();
	runner.number = application.number;
	runner.childNodes = childNodes;

	// START RUN

	//console.dir(this.parentWidget.parentWidget.runStatus.startRun);
	this.parentWidget.parentWidget.runStatus.startRun(runner);

/*
	// SET TIMED SEQUENCE: runWorkflow NEEDS A DELAY TO COMPLETE
	// EXECUTION (ON WINDOWS) BEFORE runStatus IS CALLED
	function returnWhenDone() {
		//setTimeout(continueSequence, 1000);
		return false;
	}
	function showMessage(msg) {
	}
	var seq = [
		//{ func: [showMessage, window, "repeat 10 times and pause 100ms after"], repeat: 2, pauseAfter: 100 },
			{ func: dojo.hitch(this.parentWidget, "runWorkflow"), pauseAfter: 8000 },
			{ func: [showMessage, window, "after 4000ms pause this should be seen"], pauseAfter: 1000 },
			{ func: dojo.hitch(this.parentWidget, "runStatus"), pauseAfter: 10000, repeat:100 },
			{ func: returnWhenDone } // no array, just a function to call 
		//	{ func: [dojo.hitch(//console,"log","woot"), window, "woot2?"], pauseAfter: 1000, repeat:20 }
		];

		var sequenceObject = null,
		sequenceObject = new dojox.timing.Sequence({});
*/
	}


}); // plugins.workflow.StagesMenu

