dojo.provide("plugins.workflow.Parameters");

// ALLOW THE USER TO ADD, REMOVE AND MODIFY PARAMETERS

dojo.require("dijit.TitlePane");

//dojo.require("dijit.dijit"); // optimize: load dijit layer
dojo.require("dojo.parser");
dojo.require("plugins.core.Common");

// HAS A
dojo.require("plugins.workflow.ParameterRow");

dojo.declare(
    "plugins.workflow.Parameters",
	[ dijit._Widget, dijit._Templated, plugins.core.Common ],
{

//Path to the template of this widget. 
templatePath: dojo.moduleUrl("plugins", "workflow/templates/parameters.html"),

// Calls dijit._Templated.widgetsInTemplate
widgetsInTemplate : true,

// OR USE @import IN HTML TEMPLATE
cssFiles : [ dojo.moduleUrl("plugins") + "/workflow/css/parameters.css" ],

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


startup : function () {
// DO inherited, LOAD ARGUMENTS AND ATTACH THE MAIN TAB TO THE ATTACH NODE

	// COMPLETE CONSTRUCTION OF OBJECT
	this.inherited(arguments);	 



	// ADD TO TAB CONTAINER		
	this.attachNode.addChild(this.mainTab);
	this.attachNode.selectChild(this.mainTab);
},


clear : function () {
// CLEAR THE INPUT AND OUTPUT PANES
	this.mainTab.style.visibility = "hidden";	

	this.appNameNode.innerHTML = '';

	while ( this["inputRows"].firstChild )
	{
		this["inputRows"].removeChild(this["inputRows"].firstChild);
	}
	while ( this["outputRows"].firstChild )
	{
		this["outputRows"].removeChild(this["outputRows"].firstChild);
	}
},


load : function (node, shared) {
// LOAD APPLICATION DATA INTO INPUT AND OUTPUT TITLE PANES

	// SET this.shared
	this.shared = shared;

	// SET this.stageRow
	if ( node.stageRow ) this.stageRow = node.stageRow;
	else this.stageRow = node.parentWidget;
	if ( this.stageRow == null && node.childNodes )
	{
		this.stageRow = dijit.byNode(node.childNodes[0]);
	}

	// SET this.application
	this.application = node.application;
	if ( node.application == null )
	{
		//this.application = this.parentWidget.getApplication();
		this.application = node.parentWidget.application;
	}
	//console.dir(this);
//
	// INITIALISE this.childWidgets
	if ( this.childWidgets == null ) this.childWidgets = new Array;

	// DESTROY ANY EXISTING ParameterRow CHILD WIDGETS
	while ( this.childWidgets.length > 0 )
	{
		var widget = this.childWidgets.splice(0,1)[0];
		if ( widget.destroy )	widget.destroy();
	}

	// SET PANE TITLE
	this.appNameNode.innerHTML = this.application.number + ". " + this.application.name;

	// LOAD INPUT TITLE PANE
	this.loadTitlePane("input", shared);

	// LOAD OUTPUT TITLE PANE
	this.loadTitlePane("output", shared);

	// SELECT THIS TAB PANE
	this.attachNode.selectChild(this.mainTab);

	// CALL StageRow.checkValidParameters TO CHECK THAT ALL
	// REQUIRED PARAMETER INPUTS ARE SATISFIED
	var stageRow = node.parentWidget;
	//console.dir(node);

	// DON'T force IGNORE STORED Agua.getParameterValidity DATA
	var force = false;
	stageRow.checkValidParameters(force);

	// USE THE UPDATED Agua.getParameterValidity DATA TO SET CSS 
	// CLASSES OF PARAMETER ROWS
	this.setParameterRowStyles();

},


setParameterRowStyles : function () {

	var parameterRows = this.childWidgets;
	var parameterHash = new Object;
	for ( var i = 0; i < parameterRows.length; i++ )
	{
		if ( parameterRows[i].paramtype == "input" ) 
			parameterHash[parameterRows[i].name] = parameterRows[i];
	}	
	//console.dir(parameterHash);
	//for ( var key in parameterHash )
	//{
	//}

	var stageParameters = Agua.getStageParameters(this.application);
	for ( var i = 0; i < stageParameters.length; i++ )
	{
		if ( stageParameters[i].paramtype != "input" ) continue;

		var parameterRow = parameterHash[stageParameters[i].name];

		var isValid = Agua.getParameterValidity(stageParameters[i]);
		if ( isValid == true || isValid == null )
		{
			parameterRow.setValid(parameterRow.domNode);
		}
		else
		{
			parameterRow.setInvalid(parameterRow.domNode);
		}
	}	
},


loadTitlePane : function (paneType) {
	var paneRows = paneType + "Rows";

	// CLEAR PANE
	while ( this[paneRows].firstChild )
		this[paneRows].removeChild(this[paneRows].firstChild);

	var stageObject = {
		username: this.application.username,
		project: this.application.project,
		workflow: this.application.workflow,
		name: this.application.name,
		number: this.application.number   // NB: SWITCH FROM number TO appnumber
	};

	// GET OUTPUTS FROM Agua.stageparameters	
	var parameters;
	if ( this.shared == true )
	{
		parameters = Agua.getSharedStageParameters(stageObject);
	}
	else
	{
		parameters = Agua.getStageParameters(stageObject);
	}
	parameters = this.filterByKeyValues(parameters, ["paramtype"], [paneType]);
	parameters = this.sortHasharrayByKeys(parameters, ["ordinal","name"]);

	if ( parameters == null )
	{
		return;
	}

	for ( var i = 0; i < parameters.length; i++ )
	{

		// SET parameter KEY:VALUE PAIRS
		// NB: CONVERTS ON SERVER TO HTML-SAFE (E.G., paramFunction)
		var parameter = new Object;
		for ( var key in parameters[i] )
		{
			parameter[key] = parameters[i][key];
		}

		// CONVERT PROJECT AND WORKFLOW VALUES
		if ( parameter.value == null )	parameter.value = '';
		parameter.value = parameter.value.replace(/%project%/, parameter.project);
		parameter.value = parameter.value.replace(/%workflow%/, parameter.workflow);
		parameter.value = parameter.value.replace(/%username%/, parameter.username);

		// ADD CORE LIST
		parameter.core = this.core;

		// INSTANTIATE plugins.workflow.ParameterRow
		var ParameterRow = new plugins.workflow.ParameterRow(parameter);
		this[paneRows].appendChild(ParameterRow.domNode);
		// PUSH ONTO ARRAY OF CHILD WIDGETS
		this.childWidgets.push(ParameterRow);
	}


},


checkValidInputs : function () {
// 1. CHECK VALIDITY OF ALL PARAMETERS, STORE AS this.isValid
// 2. CHANGE StageRow Style ACCORDINGLY	SET stageRow.isValid
// 3. stageRow CALLS Stages.updateValidity AND TOGGLES 
// 		runWorkflow BUTTON

	this.isValid = true;
	for ( var i = 0; i < this.childWidgets.length; i++ )
	{
		if ( this.childWidgets[i].paramtype != "input" )	continue;

		if ( this.childWidgets[i].validInput == false )
		{
			this.isValid = false;

			break;
		}
	}	

	// CALL StageRow.checkValidParameters TO CHECK THAT ALL
	// REQUIRED PARAMETER INPUTS ARE SATISFIED
	var stageRow = this.stageRow;
	if ( this.stageRow == null )	return;

	if ( this.isValid == true ) this.stageRow.setValid();
	else this.stageRow.setInvalid();


	//// CALL Stages TO UPDATE STAGE NODES AND RUN WORKFLOW BUTTON
	//// AND SET Stages.validStages FOR ENABLING RUN
	//var stages = this.parentWidget.stages;
	//
	//if ( stages != null )
	//{
	//}
}




}); // plugins.workflow.Parameters

