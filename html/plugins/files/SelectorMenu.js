
dojo.provide("plugins.files.SelectorMenu");

// WIDGET PARSER
dojo.require("dojo.parser");

// INHERITS
dojo.require("plugins.files.FileMenu");

// HAS A
dojo.require("dijit.Menu");
dojo.require("plugins.core.InputDialog");
dojo.require("plugins.core.ConfirmDialog");

dojo.declare(
    "plugins.files.SelectorMenu",
	[ plugins.files.FileMenu ],
{
	//Path to the template of this widget. 
	templatePath: dojo.moduleUrl("plugins", "files/templates/selectormenu.html"),

	//// Calls dijit._Templated.widgetsInTemplate
	//widgetsInTemplate : true,
	//

	// OR USE @import IN HTML TEMPLATE
	cssFiles : [ dojo.moduleUrl("plugins.files") + "/css/selectormenu.css" ],

	constructor : function()
	{
		//// GET INFO FROM ARGS
		//this.parentWidget = args.parentWidget;

        // LOAD SORIA AND FILEPICKER CSS
        this.loadCSS();		
	},

	postCreate : function()
	{

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


	startup : function ()
	{

		// COMPLETE CONSTRUCTION OF OBJECT
		this.inherited(arguments);	 

		// CONNECT SHORTKEYS FOR MENU
		this.setMenu();	
	},

	// ADD A NEW PROJECT USING A DIALOG BOX FOR PROJECT NAME INPUT
	select : function ()
	{

		// GET PROJECT WIDGET		
		var path = this.getPath();
		if ( path == null )
		{
			return;
		}


	},


	// ADD A NEW WORKFLOW USING A DIALOG BOX FOR WORKFLOW NAME INPUT
	add : function ()
	{


		// GET PROJECT WIDGET		
		var path = this.getPath();
		if ( path == null )
		{
			return;
		}
	}

}); // plugins.files.SelectorMenu
