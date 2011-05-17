
dojo.provide("plugins.files.FileSelectorMenu");

// WIDGET PARSER
dojo.require("dojo.parser");

// INHERITS
dojo.require("plugins.files.FileMenu");

// HAS A
dojo.require("dijit.Menu");
dojo.require("plugins.core.InputDialog");
dojo.require("plugins.core.ConfirmDialog");

dojo.declare(
    "plugins.files.FileSelectorMenu",
	[ plugins.files.FileMenu ],
{
	////}}

//Path to the template of this widget. 
templatePath: dojo.moduleUrl("plugins", "files/templates/fileselectormenu.html"),

//// Calls dijit._Templated.widgetsInTemplate
//widgetsInTemplate : true,

// OR USE @import IN HTML TEMPLATE
cssFiles : [ dojo.moduleUrl("plugins.files") + "/css/fileselectormenu.css" ],

constructor : function() {
	//// GET INFO FROM ARGS
	//this.parentWidget = args.parentWidget;

	// LOAD SORIA AND FILEPICKER CSS
	this.loadCSS();		
},

postCreate : function() {

	//// SET INPUT DIALOG
	//this.setInputDialog();
	//
	//// SET CONFIRM DIALOG
	//this.setConfirmDialog();

	// SET LABEL
	this.setTitle("Selector Menu");

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

select : function (event) {
// ADD A NEW PROJECT USING A DIALOG BOX FOR PROJECT NAME INPUT
	//////console.dir(event.target);

	// GET PROJECT WIDGET		
	var location = this.getPath();
	if ( location == null )
	{
		return;
	}

	var filename = this.menu.currentTarget.innerHTML;

	this.selectCallback(filename, location, this.type, this.parameterWidget);
	this.hide();
},


add : function () {
// ADD A NEW WORKFLOW USING A DIALOG BOX FOR WORKFLOW NAME INPUT


	// GET PROJECT WIDGET		
	var location = this.getPath();
	if ( location == null )
	{
		return;
	}

	var filename = this.menu.currentTarget.innerHTML;

	this.addCallback(filename, location, this.type, this.parameterWidget);
	this.hide();
}

}); // plugins.files.FileSelectorMenu
