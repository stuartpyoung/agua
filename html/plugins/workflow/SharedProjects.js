dojo.provide("plugins.workflow.SharedProjects");

/* CLASS SUMMARY: DISPLAY SHARED PROJECTS AND WORKFLOWS

   ALLOW THE USER TO:

		-	VIEW THE PARAMETERS OF SHARED WORKFLOWS 

		-	VIEW THE OUTPUT FILES OF SHARED WORKFLOWS

		-	COPY SHARED WORKFLOWS TO USER'S OWN PROJECTS

   TO DO:

	copyProject
	sharedStatus

*/

dojo.require("dijit.layout.ContentPane");
dojo.require("dojo.parser");
dojo.require("plugins.core.Common");
dojo.require("plugins.workflow.Stages");

dojo.declare(
    "plugins.workflow.SharedProjects",
	[ plugins.workflow.Stages ],
{
//Path to the template of this widget. 
templatePath: dojo.moduleUrl("plugins", "workflow/templates/sharedprojects.html"),

// Calls dijit._Templated.widgetsInTemplate
widgetsInTemplate : true,

// OR USE @import IN HTML TEMPLATE
// **** LOADED IN Workflow.js ****
//cssFiles : [ dojo.moduleUrl("plugins") + "/workflow/css/shared.css" ],

shared : true,

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

	//// COMPLETE CONSTRUCTION OF OBJECT
	//this.inherited(arguments);	 


	// ADD TO TAB CONTAINER		
	this.tabContainer.addChild(this.mainTab);
	this.tabContainer.selectChild(this.mainTab);

	// SET SELECTIVE DIALOG FOR copyWorkflow	
	this.setSelectiveDialog();

	// SET INTERACTIVE DIALOG FOR copyProject
	this.setInteractiveDialog();

	// CREATE SOURCE MENU
	this.setContextMenu();

	// SET SHARING USER NAMES COMBO
	this.setUsernameCombo();
},


enableRunButton : function () {
},

disableRunButton : function () {
},

setUsernameCombo : function () {

	this.inherited(arguments);

	// DOJO.CONNECT TO CHANGE THE workflowCombo
	var thisObject = this;
	dojo.connect(this.usernameCombo, "onchange", dojo.hitch(function(event) {
		thisObject.setProjectCombo(event);
	}));

	var username = this.usernameCombo.get('value');

	// SET THE PROJECT COMBO
	this.setProjectCombo(username);
},


setProjectCombo : function (username, projectName, workflowName) {

	this.setSharedProjectCombo(username, projectName, workflowName);

	// DOJO.CONNECT TO CHANGE THE workflowCombo
	var thisObject = this;
	dojo.connect(this.projectCombo, "onchange", dojo.hitch(function(event) {

		var username = thisObject.usernameCombo.get('value');
		var projectName = thisObject.projectCombo.get('value');

		thisObject.setWorkflowCombo(username, projectName);
	}));

	if ( projectName == null || ! projectName )
	{
		projectName = this.projectCombo.get('value');
	}

	// SET THE PROJECT COMBO
	this.setWorkflowCombo(username, projectName);
},


// SET THE PROJECT COMBO
setWorkflowCombo : function (username, projectName, workflowName) {

	this.setSharedWorkflowCombo(username, projectName, workflowName);

	// DOJO.CONNECT TO CHANGE THE workflowCombo
	var thisObject = this;
	dojo.connect(this.workflowCombo, "onchange", dojo.hitch(function(event) {

		var username = thisObject.usernameCombo.get('value');
		var projectName = thisObject.projectCombo.get('value');
		var workflowName = thisObject.workflowCombo.get('value');
		thisObject.setDropTarget(username, projectName, workflowName);
	}));


	if ( projectName == null || ! projectName )
		projectName = this.projectCombo.get('value');

	if ( workflowName == null || ! workflowName )
		workflowName = this.workflowCombo.get('value');

	this.setDropTarget(username, projectName, workflowName);
},



setDropTarget : function (username, projectName, workflowName) {

	// SET THE DROP TARGET
	var sharedStages  = Agua.getSharedStagesByWorkflow(username, projectName, workflowName);
	if ( sharedStages == null )
	{
		return;
	}

	this._setDropTarget(sharedStages);
},


assumeFocus : function () {

	var username = thisObject.usernameCombo.get('value');
	var projectName = thisObject.projectCombo.get('value');
	var workflowName = thisObject.workflowCombo.get('value');
	this.setDropTarget(username, projectName, workflowName);
},


// LOAD DATA INTO INFO PANE FROM THE APPLICATION ASSOCIATED WITH THIS NODE
// OVERLOAD THIS TO PASS ADDITIONAL ARGUMENTS TO Parameters.load()
loadParametersPane : function (node) {

	// WARN AND QUIT IF NO NODE PASSED, E.G., IF WORKFLOW HAS NO STAGES
	if ( node == null )
	{
		return;
	}

	if ( this.parentWidget.parameterPane != null )
		this.parentWidget.parameterPane.load(node, this.shared);
},

_copyWorkflow : function (sourceProject, sourceWorkflow, targetProject, targetWorkflow, copyFiles) {

	var sourceUser = this.usernameCombo.get('value');
	var targetUser = Agua.cookie('username');


	// ADD PROJECT
	Agua.copyWorkflow(sourceUser, sourceProject, sourceWorkflow, targetUser, targetProject, targetWorkflow, copyFiles);

},


_copyProject : function (sourceProject, targetProject, copyFiles) {

	var targetUser = Agua.cookie('username');
	var sourceUser =  thisObject.usernameCombo.get('value');

	// ADD PROJECT
	Agua.copyProject(sourceUser, sourceProject, targetUser, targetProject, copyFiles);
}




}); // plugins.workflow.SharedProjects

/*
_copyWorkflow : function (sourceProject, sourceWorkflow, targetProject, targetWorkflow, copyFiles) {

		var targetUser = Agua.cookie('username');
		var sourceUser =  thisObject.usernameCombo.get('value');

		// ADD PROJECT
		Agua.copyWorkflow(sourceUser, sourceProject, sourceWorkflow, targetUser, targetProject, targetWorkflow, copyFiles);
},
*/
