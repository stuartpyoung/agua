
dojo.provide("plugins.project.Controller");

// OBJECT:  plugins.project.Controller
// PURPOSE: GENERATE AND MANAGE Project PANES

// INHERITS
dojo.require("plugins.core.Common");

// GLOBAL ADMIN CONTROLLER VARIABLE
var projectController;

// HAS
dojo.require("plugins.project.Project");

dojo.declare( "plugins.project.Controller",
	[ dijit._Widget, dijit._Templated, plugins.core.Common ],
{
	// PANE ID 
	paneId : null,

	//Path to the template of this widget. 
	templatePath: dojo.moduleUrl("plugins", "project/templates/controller.html"),

	// Calls dijit._Templated.widgetsInTemplate
	widgetsInTemplate : true,

	// CSS FILES
	cssFiles : ["plugins/project/css/controller.css"],

	// CONSTRUCTOR	
	constructor : function(args) {

		// LOAD CSS FOR BUTTON
		this.loadCSS();
	},


	postCreate : function()
	{

		this.startup();
	},


	startup : function ()
	{

		this.inherited(arguments);

		// ADD ADMIN TAB TO TAB CONTAINER		
		Agua.toolbar.addChild(this.projectButton);

		// SET ADMIN BUTTON LISTENER
		var listener = dojo.connect(this.projectButton, "onClick", this, "createTab");

		// CREATE TAB
		this.createTab();		
	},


	createTab : function ()
	{

		var project = new plugins.project.Project( { paneId: "projectPane" } );
		this.paneId++;
	},


}); // end of Controller
