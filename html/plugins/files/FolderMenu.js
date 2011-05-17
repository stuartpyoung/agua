
dojo.provide("plugins.files.FolderMenu");

// WIDGET PARSER
dojo.require("dojo.parser");

// INHERITS
dojo.require("plugins.files.FileMenu");

// HAS A
dojo.require("dijit.Menu");

dojo.declare(
    "plugins.files.FolderMenu",
	[ plugins.files.FileMenu ],
{

	constructor : function(args)
	{
		// GET INFO FROM ARGS
		//this.parentWidget = args.parentWidget;

        // LOAD SORIA AND FILEPICKER CSS
        this.loadCSS();		
	},

	postCreate : function()
	{

		// DISABLE DOWNLOAD
		this.downloadNode.destroy();

		// SET LABEL
		this.setTitle("Folder Options");

		// SET INPUT DIALOG
		this.setInputDialog();

		// SET CONFIRM DIALOG
		this.setConfirmDialog();

		// CONNECT SHORTKEYS FOR MENU
		this.setMenu();

		// SET THE UPLOAD OBJECT
		this.setUploader();

		// DO STARTUP
		this.startup();
	}



}); // plugins.files.FolderMenu
