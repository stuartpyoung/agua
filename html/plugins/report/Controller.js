
dojo.provide("plugins.report.Controller");

// OBJECT:  plugins.report.Controller
// PURPOSE: GENERATE AND MANAGE Report PANES

// INHERITS
//dojo.require("plugins.core.WidgetFramework");

// GLOBAL ADMIN CONTROLLER VARIABLE
var reportController;

// INHERITS
dojo.require("plugins.core.Common");

// HAS
dojo.require("plugins.report.Report");


dojo.declare( "plugins.report.Controller",
	[ dijit._Widget, dijit._Templated, plugins.core.Common ],
{
	//Path to the template of this widget. 
	templatePath: dojo.moduleUrl("plugins", "report/templates/controller.html"),

	// Calls dijit._Templated.widgetsInTemplate
	widgetsInTemplate : true,

	// CSS FILE FOR BUTTON STYLING
	cssFiles : [ "plugins/report/css/controller.css" ],

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
		this.reportButton = Agua.addToolbarButton("Report");
		Agua.toolbar.addChild(this.reportButton);

		// SET ADMIN BUTTON LISTENER
		var listener = dojo.connect(this.reportButton, "onClick", this, "createTab");

		// CREATE TAB
		//this.createTab();		

	},


	createTab : function (args)
	{

		if ( args == null ) args = new Object;
		args.attachWidget = Agua.tabs;
		this.tabPanes.push(new plugins.report.Report(args));
	}

}); // end of Controller


//dojo.addOnLoad(
//	function()
//	{
//		if ( Agua.loginController == null )
//		{
//			
//			Agua.reportController = new plugins.report.Controller();
//			
//			// DEBUG
//			Agua.reportController.createTab();
//		}	
//	}
//);

