
dojo.provide("plugins.files.WorkflowMenu");

// WIDGET PARSER
dojo.require("dojo.parser");

// INHERITS
dojo.require("plugins.files.FileMenu");

// HAS A
dojo.require("dijit.Menu");
//dojo.require("plugins.core.InputDialog");
//dojo.require("plugins.core.InteractiveDialog");
//dojo.require("plugins.core.ConfirmDialog");

dojo.declare(
    "plugins.files.WorkflowMenu",
	[ plugins.files.FileMenu ],
{
	/////}}	

//Path to the template of this widget. 
templatePath: dojo.moduleUrl("plugins", "files/templates/workflowmenu.html"),

//// Calls dijit._Templated.widgetsInTemplate
//widgetsInTemplate : true,
//

// OR USE @import IN HTML TEMPLATE
cssFiles : [ dojo.moduleUrl("plugins.files") + "/css/workflowmenu.css" ],

constructor : function() {

	//// GET INFO FROM ARGS
	//this.parentWidget = args.parentWidget;

	// LOAD SORIA AND FILEPICKER CSS
	this.loadCSS();		
},

postCreate : function() {

	// SET INPUT DIALOG
	this.setInputDialog();

	// SET INTERACTIVE DIALOG
	this.setInteractiveDialog();

	// SET CONFIRM DIALOG
	this.setConfirmDialog();

	// SET LABEL
	this.setTitle("Workflow Menu");

	// CONNECT SHORTKEYS FOR MENU
	this.setMenu();

	// DO STARTUP
	this.startup();
},


startup : function () {

	// COMPLETE CONSTRUCTION OF OBJECT
	this.inherited(arguments);	 

	// CONNECT SHORTKEYS FOR MENU
	this.setMenu();	
},


setTitle : function (title) {
// NO TITLE - DO NOTHING

//this.titleNode.containerNode.innerHTML = title;
},

copyWorkflow : function () {

},

renameWorkflow : function () {

	// GET PROJECT WIDGET		
	var projectWidget = this.getProjectWidget();
	if ( projectWidget == null )	return;

	var projectName = this.getProjectName();

	var originalWorkflowName = this.getWorkflowName();

	// SET TITLE AND MESSAGE
	var title = "Rename workflow '" + originalWorkflowName + "'";
	var message = "Please enter new name";

	// CALLBACKS
	var cancelCallback = function (){
	};
	var enterCallback = dojo.hitch(this, function (newWorkflowName)
		{

			// SANITY CHECK
			if ( newWorkflowName == null )	return;
			newWorkflowName = newWorkflowName.replace(/\s+/, '');
			if ( newWorkflowName == '' )	return;

			// CHECK IF NAME EXISTS ALREADY
			if ( Agua.isWorkflow(projectName, newWorkflowName) == true )
			{
				return;
			}


			// GET ORIGINAL WORKFLOW OBJECT
			var workflowObject = Agua.getWorkflow(projectName, originalWorkflowName);

			// QUIT IF WORKFLOW DOES NOT EXIST
			if ( workflowObject == null )
			{
				return;
			}

			// ADD PROJECT
			Agua.renameWorkflow(workflowObject, newWorkflowName);

			// RELOAD THE PROJECTS TAB
			setTimeout(function(thisObj) { projectWidget.reload(); }, 1000, this);
		}
	);		

	// SHOW THE DIALOG
	this.loadInputDialog(title, message, enterCallback, cancelCallback);

},

openWorkflow : function () {

},

newProject : function () {
// ADD A NEW PROJECT USING A DIALOG BOX FOR PROJECT NAME INPUT

	// GET PROJECT WIDGET		
	var projectWidget = this.getProjectWidget();
	if ( projectWidget == null )
	{
		return;
	}


	// CALLBACKS
	var cancelCallback = function (){
	};
	var enterCallback = dojo.hitch(this, function (projectName)
		{

			// SANITY CHECK
			if ( projectName == null )	return;
			if ( projectName == '' )	return;

			projectName = projectName.replace(/\s+/, '');

			// THEN ADD NEW PROJECT OBJECT TO Agua.projects ARRAY
			var projectObject = new Object;
			projectObject.name = projectName;

			if ( Agua.isProject(projectName) == true )
			{
				return;
			}

			// ADD PROJECT
			Agua.addProject(projectObject);

			// RELOAD THE PROJECTS TAB
			setTimeout(function(thisObj) { projectWidget.reload(); }, 1000, this);
		}
	);		

	var title = "New Project";
	var message = "Please enter project name";
	this.loadInputDialog(title, message, enterCallback, cancelCallback);
},


newWorkflow : function () {
// ADD A NEW WORKFLOW USING A DIALOG BOX FOR WORKFLOW NAME INPUT

	// GET PROJECT WIDGET		
	var projectWidget = this.getProjectWidget();
	if ( projectWidget == null )	return;

	var projectName = this.getProjectName();

	var interactiveDialog = this.interactiveDialog;

	// SET TITLE AND MESSAGE
	var title = "New Workflow";
	var message = "Please enter workflow name";

	// CALLBACKS
	var cancelCallback = function (){
	};

	// CALLBACK FORMAT:
	// this.dialog.enterCallback(input, checked);	
	var enterCallback = dojo.hitch(this, function (workflowName, undefined )
		{

			// SANITY CHECK
			if ( workflowName == null || workflowName == '' )	return;
			workflowName = workflowName.replace(/\s+/, '');

			// QUIT IF WORKFLOW EXISTS ALREADY
			if ( Agua.isWorkflow(projectName, workflowName) == true )
			{

				interactiveDialog.messageNode.innerHTML = "Workflow already exists in '" + projectName + "'";
				return;
			}
			else {

				interactiveDialog.messageNode.innerHTML = "Creating workflow";
				interactiveDialog.close();
			}

			// ADD PROJECT
			Agua.addWorkflow(projectName, workflowName);

			// RELOAD THE PROJECTS TAB
			setTimeout(function(thisObj) { projectWidget.reload(); }, 1000, this);
		}
	);		

	// SHOW THE DIALOG
	this.loadInteractiveDialog(title, message, enterCallback, cancelCallback);
},

deleteWorkflow : function () {
// DELETE A WORKFLOW USING A DIALOG BOX FOR CONFIRMATION BY THE USER

	// GET PROJECT WIDGET		
	var projectWidget = this.getProjectWidget();
	if ( projectWidget == null )	return;

	var projectName = this.getProjectName();

	var workflowName = this.getWorkflowName();

	// CALLBACKS
	var noCallback = function (){
	};
	var yesCallback = dojo.hitch(this, function ()
		{

			// SANITY CHECK
			if ( workflowName == null || workflowName == '' )	return;
			workflowName = workflowName.replace(/\s+/, '');

			// QUIT IF WORKFLOW EXISTS ALREADY
			if ( Agua.isWorkflow(projectName, workflowName) == false )
			{
				return;
			}

			// REMOVE WORKFLOW
			var workflowObject = new Object;
			workflowObject.project = projectName;
			workflowObject.name = workflowName;
			Agua.removeWorkflow(workflowObject);

			// RELOAD THE PROJECTS TAB
			setTimeout(function(thisObj) { projectWidget.reload(); }, 1000, this);
		}
	);		

	// SET TITLE AND MESSAGE
	var title = "Delete workflow '" + workflowName + "'?";
	var message = "All files and data will be destroyed";

	// SHOW THE DIALOG
	this.loadConfirmDialog(title, message, yesCallback, noCallback);
},

deleteProject : function () {
// DELETE A PROJECT USING A DIALOG BOX FOR CONFIRMATION BY THE USER

	// GET PROJECT WIDGET		
	var projectWidget = this.getProjectWidget();
	if ( projectWidget == null )	return;

	var projectName = this.getProjectName();

	// CALLBACKS
	var noCallback = function (){
	};
	var yesCallback = dojo.hitch(this, function ()
		{

			// SANITY CHECK
			if ( Agua.isProject(projectName) == false )
			{
				return;
			}

			// REMOVE PROJECT
			Agua.removeProject({ name: projectName });

			// RELOAD THE PROJECTS TAB
			setTimeout(function(thisObj) { projectWidget.reload(); }, 1000, this);
		}
	);		

	// SET TITLE AND MESSAGE
	var title = "Delete project '" + projectName + "'?";
	var message = "All workflows and data will be destroyed";

	// SHOW THE DIALOG
	this.loadConfirmDialog(title, message, yesCallback, noCallback);

},

getProjectWidget : function () {
// RETURN THE PROJECT TAB WIDGET CONTAINING THIS FILE DRAG OBJECT

	// SANITY		
	if ( this.menu.currentTarget == null )	return null;

	// GET THE PROJECT WIDGET
	var item = this.menu.currentTarget.item;
	var widget = dijit.getEnclosingWidget(this.menu.currentTarget);
	var projectWidget = widget.parentWidget.parentWidget;

	return projectWidget;
},

getProjectName : function () {
// RETURN THE PROJECT NAME FOR THIS FILE DRAG OBJECT

	// SANITY		
	if ( this.menu.currentTarget == null )	return null;

	// GET THE PROJECT WIDGET
	var item = this.menu.currentTarget.item;
	var widget = dijit.getEnclosingWidget(this.menu.currentTarget);
	var projectName = widget.path;

	return projectName;
},

getWorkflowName : function () {
// RETURN THE WORKFLOW NAME FOR THIS GROUP DRAG PANE OBJECT

	// SANITY		
	if ( this.menu.currentTarget == null )	return null;

	// GET THE PROJECT WIDGET
	var item = this.menu.currentTarget.item;
	var workflowName = item.path;

	return workflowName;
}


}); // plugins.files.WorkflowMenu
