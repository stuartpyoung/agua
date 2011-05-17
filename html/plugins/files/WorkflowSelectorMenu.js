
dojo.provide("plugins.files.WorkflowSelectorMenu");

// WIDGET PARSER
dojo.require("dojo.parser");

// INHERITS
dojo.require("plugins.files.FileMenu");

// HAS A
dojo.require("dijit.Menu");

dojo.declare(
    "plugins.files.WorkflowSelectorMenu",
	[ plugins.files.FileMenu ],
{
	//Path to the template of this widget. 
	templatePath: dojo.moduleUrl("plugins", "files/templates/workflowselectormenu.html"),


	//// Calls dijit._Templated.widgetsInTemplate
	//widgetsInTemplate : true,
	//
	// OR USE @import IN HTML TEMPLATE
	cssFiles : [ dojo.moduleUrl("plugins.files") + "/css/workflowselectormenu.css" ],
	//
	//// PARENT WIDGET
	//parentWidget : null,

	//////}

constructor : function(args)  {
	// GET INFO FROM ARGS
	//this.parentWidget = args.parentWidget;

	// LOAD SORIA AND FILEPICKER CSS
	this.loadCSS();		
},

postCreate : function() {


	// SET LABEL
	this.setTitle("Folder Options");

	// SET INPUT DIALOG
	this.setInputDialog();

	// CONNECT SHORTKEYS FOR MENU
	this.setMenu();

	// DO STARTUP
	this.startup();
}



}); // plugins.files.WorkflowSelectorMenu
