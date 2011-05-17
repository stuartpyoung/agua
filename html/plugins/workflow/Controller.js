dojo.provide("plugins.workflow.Controller");

// OBJECT:  plugins.workflow.Controller
// PURPOSE: GENERATE AND MANAGE Workflow PANES

// INHERITS
dojo.require("plugins.core.Common");

// HAS
dojo.require("plugins.workflow.Workflow");

dojo.declare( "plugins.workflow.Controller",
	[ dijit._Widget, dijit._Templated, plugins.core.Common ],
{

	//Path to the template of this widget. 
	templatePath: dojo.moduleUrl("plugins", "workflow/templates/controller.html"),

	// Calls dijit._Templated.widgetsInTemplate
	widgetsInTemplate : true,

	// CSS FILE FOR BUTTON STYLING
	cssFiles : [ dojo.moduleUrl("plugins") + "/workflow/css/controller.css", "plugins/workflow/css/workflow.css" ],

	// ARRAY OF TAB PANES
	tabPanes : [],

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
		this.workflowButton = Agua.addToolbarButton("Workflow");
		Agua.toolbar.addChild(this.workflowButton);

		// SET ADMIN BUTTON LISTENER
		var listener = dojo.connect(this.workflowButton, "onClick", this, "createTab");

		// CREATE TAB
		this.createTab();		

	},


	createTab : function (args)
	{

		if ( args == null ) args = new Object;
		args.attachWidget = Agua.tabs;
		this.tabPanes.push(new plugins.workflow.Workflow(args));
	}


}); // end of Controller

