dojo.provide("plugins.view.Controller");

// OBJECT:  plugins.view.Controller
// PURPOSE: GENERATE AND MANAGE View PANES

// INHERITS
dojo.require("plugins.core.Common");

// GLOBAL ADMIN CONTROLLER VARIABLE
var viewController;

// HAS
dojo.require("plugins.view.View");

dojo.declare( "plugins.view.Controller",
	[ dijit._Widget, dijit._Templated, plugins.core.Common ],
{
	/////}}

//Path to the template of this widget. 
templatePath: dojo.moduleUrl("", "../plugins/view/templates/controller.html"),

// Calls dijit._Templated.widgetsInTemplate
widgetsInTemplate : true,

// CSS FILE FOR BUTTON STYLING
cssFiles : [ "plugins/view/css/controller.css" ],

// ARRAY OF TAB PANES
tabPanes : [],

// CONSTRUCTOR	
constructor : function(args) {

	// LOAD CSS FOR BUTTON
	this.loadCSS();

},


postCreate : function() {

	this.startup();
},


startup : function () {

	this.inherited(arguments);

	// ADD TO TAB CONTAINER		
	this.viewButton = Agua.addToolbarButton("View");

	// SET BUTTON LISTENER
	var listener = dojo.connect(this.viewButton, "onClick", this, "createTab");

	// AUTOMATICALLY CREATE TAB ON LAUNCH
	this.createTab();
},


createTab : function (args) {

	if ( args == null ) args = new Object;
	args.attachWidget = Agua.tabs;
	this.tabPanes.push(new plugins.view.View(args));
}


}); // end of Controller


dojo.addOnLoad(
	function()
	{

		if ( Agua.loginWindow == null )
		{
			//Agua.viewController = new plugins.view.Controller();

			// DEBUG
			//Agua.viewController.createTab();
		}
	}
);

