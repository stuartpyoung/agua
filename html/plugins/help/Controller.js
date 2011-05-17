
dojo.provide("plugins.help.Controller");

// OBJECT:  plugins.help.Controller
// PURPOSE: GENERATE AND MANAGE Help PANES

// INHERITS
dojo.require("plugins.core.Common");

// GLOBAL ADMIN CONTROLLER VARIABLE
var helpController;

// HAS
dojo.require("plugins.help.Help");

dojo.declare( "plugins.help.Controller",
	[ dijit._Widget, dijit._Templated, plugins.core.Common ],
{
	//Path to the template of this widget. 
	templatePath: dojo.moduleUrl("plugins", "help/templates/controller.html"),

	// Calls dijit._Templated.widgetsInTemplate
	widgetsInTemplate : true,

	// CSS FILE FOR BUTTON STYLING
	cssFiles : [ "plugins/help/css/controller.css" ],

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
		this.helpButton = Agua.addToolbarButton("Help");
		Agua.toolbar.addChild(this.helpButton);

		// SET ADMIN BUTTON LISTENER
		var listener = dojo.connect(this.helpButton, "onClick", this, "createTab");
	},


	createTab : function (args)
	{

		if ( args == null ) args = new Object;
		args.attachWidget = Agua.tabs;
		this.tabPanes.push(new plugins.help.Help(args));
	}


}); // end of Controller

