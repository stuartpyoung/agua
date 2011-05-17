dojo.provide("plugins.workflow.Shared");

// DISPLAY PROJECTS/WORKFLOWS THAT HAVE BEEN SHARED WITH THE USER

dojo.require("dijit.layout.ContentPane");

//dojo.require("dijit.dijit"); // optimize: load dijit layer
dojo.require("dojo.parser");
dojo.require("plugins.core.Common");

// HAS A
//dojo.require("plugins.workflow.SharedRow");

dojo.declare(
    "plugins.workflow.Shared",
	[ dijit._Widget, dijit._Templated, plugins.core.Common ],
{
//Path to the template of this widget. 
templatePath: dojo.moduleUrl("plugins", "workflow/templates/shared.html"),

// Calls dijit._Templated.widgetsInTemplate
widgetsInTemplate : true,

// OR USE @import IN HTML TEMPLATE
// **** LOADED IN Workflow.js ****
//cssFiles : [ "plugins/workflow/css/shared.css" ],

// PARENT WIDGET
parentWidget : null,

// ARRAY OF CHILD WIDGETS
childWidgets : null,

// isVALID BOOLEAN: ALL PARAMETERS ARE VALID
isValid : null,

// CORE WORKFLOW OBJECTS
core : null,

/////}

constructor : function(args) {
	this.core = args.core;

	// LOAD CSS
	this.loadCSS();		
},

postCreate : function() {

	this.startup();
},


// DO inherited, LOAD ARGUMENTS AND ATTACH THE MAIN TAB TO THE ATTACH NODE
startup : function () {

	// COMPLETE CONSTRUCTION OF OBJECT
	this.inherited(arguments);	 


	// ADD TO TAB CONTAINER		
	this.attachNode.addChild(this.mainTab);
	this.attachNode.selectChild(this.mainTab);

	// SET SHARING USER NAMES COMBO
	this.setUsernameCombo();
},

setUsernameCombo : function () {

	this.inherited(arguments);

	// DOJO.CONNECT TO CHANGE THE workflowCombo
	var thisObject = this;
	dojo.connect(this.usernameCombo, "onChange", function(event) {
		thisObject.setSharedProjectCombo(event);
	});

	var username = this.usernameCombo.getValue();

	// SET THE PROJECT COMBO
	this.setSharedProjectCombo(username);
},


setSharedProjectCombo : function (username, projectName, workflowName) {

	this.inherited(arguments);

	// DOJO.CONNECT TO CHANGE THE workflowCombo
	var thisObject = this;
	dojo.connect(this.sharedProjectCombo, "onChange", function(event) {

		var userName = thisObject.usernameCombo.getValue();
		var projectName = thisObject.sharedProjectCombo.getValue();

		thisObject.setSharedWorkflowCombo(userName, projectName);
	});

	if ( projectName == null || ! projectName )
	{
		projectName = this.sharedProjectCombo.getValue();
	}

	// SET THE PROJECT COMBO
	this.setSharedWorkflowCombo(username, projectName);
},


// SET THE PROJECT COMBO
setSharedWorkflowCombo : function (username, projectName, workflowName) {

	this.inherited(arguments);

	// DOJO.CONNECT TO CHANGE THE workflowCombo
	var thisObject = this;
	dojo.connect(this.sharedWorkflowCombo, "onChange", function(event) {

		var userName = thisObject.usernameCombo.getValue();
		var projectName = thisObject.sharedProjectCombo.getValue();
		var workflowName = thisObject.sharedWorkflowCombo.getValue();
		var sharedStages  = Agua.getSharedStagesByWorkflow(username, projectName, workflowName);
		thisObject.setSharedDropTarget(sharedStages);
	});


	if ( projectName == null || ! projectName )
		projectName = this.sharedProjectCombo.getValue();

	if ( workflowName == null || ! workflowName )
		workflowName = this.sharedWorkflowCombo.getValue();

	// SET THE DROP TARGET
	var sharedStages  = Agua.getSharedStagesByWorkflow(username, projectName, workflowName);
	if ( sharedStages == null )
	{
		return;
	}

	this.setSharedDropTarget(sharedStages);
},




setSharedDropTarget : function (stages) {

	// GET APPLICATIONS FOR THIS WORKFLOW IN THIS PROJECT
	if ( stages == null || ! stages )
	{
		return;
	}

	// GENERATE LIST OF TARGET APPLICATIONS FOR DRAG AND DROP
	// FROM WORKFLOW APPLICATIONS LIST
	var dataArray = new Array;
	for ( var i = 0; i < stages.length; i++ )
	{
		var application = stages[i];


		var applicationName = application.name;
		if ( applicationName )
		{
			var hash = new Object;
			hash.data = applicationName;
			hash.type = [ "isDraggable" ];
			dataArray.push(hash);
		}
	}

	this.sharedDropTarget = new dojo.dnd.Target( this.sharedDropTargetContainer, { accept: ["isDraggable"] } );
	dojo.addClass(this.sharedDropTarget, 'dropTarget');

	// DELETE ALL NODES	(NEEDED TO CLEAN UP REMAINING NODES FROM PREVIOUS WORKFLOW)	
	var selectedNodes = this.sharedDropTarget.selectAll();
	this.sharedDropTarget.deleteSelectedNodes(selectedNodes);

	this.sharedDropTarget.insertNodes(false, dataArray);


	// USE IN dojo.connect AND OVERRIDE OF onDndDrop BELOW
	var thisObject = this;


	// SET NODE CHARACTERISTICS - ONCLICK, CLASS, ETC.
	allNodes = this.sharedDropTarget.getAllNodes();
	for ( var i = 0; i < allNodes.length; i++ )
	{

		var node = allNodes[i];
		var nodeClass = dataArray[i].type;
		var applicationName = node.innerHTML;

		// ADD TARGET MENU
		//this.targetMenu.bindDomNode(node);

		// ADD workflow NAMESPACE CLASS TO NODE
		dojo.addClass(node, 'workflow');

		// ADD CLASS FROM type TO NODE
		dojo.addClass(node, nodeClass);

		// ADD infoId TO NODE
		node.setAttribute('infoId', thisObject.infoId);

		// ADD application TO NODE
		node.application = stages[i];

		// ADD application TO NODE
		node.application = stages[i];

		// NO DESCRIPTION OR NOTES FOR NOW IN STAGE
		stages[i].description = '';
		stages[i].notes = '';

		// INSTANTIATE SOURCE ROW 
		var appRow = new plugins.workflow.StageRow(stages[i]);

		// CLEAR NODE CONTENT
		node.innerHTML = '';

		// APPEND TO NODE
		node.appendChild(stageRow.domNode);

		// ADD stageRow AS node.stagerow ATTRIBUTE FOR ACCESS LATER
		// WHEN CALLING Workflow.loadParametersPane SO THAT THE CORRECT
		// StageRow HAS ITS validInputs SET ACCORDING TO THE OUTCOME
		// OF Workflow.loadParametersPane
		node.stageRow = stageRow;

		// HACK:
		//
		// SET appRow AS node.parentWidget FOR LATER RESETTING OF
		// number ON REMOVAL OR INSERTION OF NODES
		//
		// REM: remove ONCLICK BUBBLES ON appRow.name NODE RATHER THAN ON node. 
		// I.E., CONTRARY TO DESIRED, this.name IS THE TARGET INSTEAD OF THE node.
		//
		// ALSO ADDED this.name.parentWidget = this IN AppRow.startup()
		node.parentWidget = appRow;

		// SHOW APPLICATION INFO WHEN CLICKED
		dojo.connect(node, "onclick",  null, function(event)
			{
				var shared = true;

				thisObject.parentWidget.parentWidget.loadParametersPane(event.target, shared);
			}
		);
	}	// END OF for ( var i = 0; i < allNodes.length; i++ )

	// DO loadParametersPane REGARDLESS OF WHETHER OR NOT THERE ARE ANY NODES
	// (IT WILL CLEAR THE INFO PANE IF THERE ARE NO NODES)
	if ( allNodes.length != null )
	{
		if ( allNodes.length == 0 )
		{
			this.parentWidget.parentWidget.loadParametersPane();
		}
		else
		{

			var shared = true;
			this.parentWidget.parentWidget.loadParametersPane(allNodes[0], shared);
		}
	}
}



}); // plugins.workflow.Shared
