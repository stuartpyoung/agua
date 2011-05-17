dojo.provide( "plugins.core.ConfirmDialog");

// HAS A
dojo.require("dijit.Dialog");
dojo.require("dijit.form.Button");

// INHERITS
dojo.require("plugins.core.Common");


dojo.declare( "plugins.core.ConfirmDialog",
	[ dijit._Widget, dijit._Templated, plugins.core.Common ],
{
	//Path to the template of this widget. 
	templatePath: dojo.moduleUrl("plugins", "core/templates/confirmdialog.html"),

	// OR USE @import IN HTML TEMPLATE
	cssFiles : [ dojo.moduleUrl("plugins.core") + "/css/confirmdialog.css" ],

	// Calls dijit._Templated.widgetsInTemplate
	widgetsInTemplate : true,

	// PARENT plugins.workflow.Apps WIDGET
	parentWidget : null,

	// APPLICATION OBJECT
	application : null,

	// DIALOG TITLE
	title: null,

	// DISPLAYED MESSAGE 
	message : null,

	constructor : function(args)
	{

		this.title 				=	args.title;
		this.message 			=	args.message;
		this.parentWidget 		=	args.parentWidget;
		this.yesCallback 		=	args.yesCallback;
		this.noCallback 		=	args.noCallback;

		// LOAD CSS
        this.loadCSS();
	},

	getApplication : function ()
	{
		return this.application;
	},

	postCreate : function()
	{

		this.startup();
	},

	startup : function ()
	{

		this.inherited(arguments);

		// SET UP DIALOG
		this.setDialogue();

		// ADD CSS NAMESPACE CLASS
		dojo.addClass(this.dialog.containerNode, "confirmDialog");
		dojo.addClass(this.dialog.titleNode, "confirmDialog");
		dojo.addClass(this.dialog.closeButtonNode, "confirmDialog");

		// REMOVE CLOSE BUTTON NODE
		this.dialog.closeButtonNode.setAttribute('display', 'none');
		this.dialog.closeButtonNode.setAttribute('visibility', 'hidden');
	},

	// SHOW THE DIALOGUE
	show: function ()
	{
		this.dialog.show();
	},

	// HIDE THE DIALOGUE
	hide: function ()
	{
		this.dialog.hide();
	},

	doYes : function(type)
	{

		// DO CALLBACK
		this.dialog.yesCallback();

		// HIDE
		this.dialog.hide();
	},

	doNo : function()
	{

		// DO CALLBACK
		this.dialog.noCallback();

		// HIDE
		this.dialog.hide();
	},


	// LOAD THE DIALOGUE VALUES
	load : function (args)
	{

		if ( args.title == null )	{	args.title = "";	}
		this.dialog.titleNode.innerHTML	=	args.title;
		this.messageNode.innerHTML	=	args.message;
		this.dialog.yesCallback		=	args.yesCallback;
		this.dialog.noCallback		=	args.noCallback


		this.show();
	},


	// APPEND TO DOCUMENT BODY
	setDialogue : function () {

		// APPEND DIALOG TO DOCUMENT
		//this.dialog.title = title;
		document.body.appendChild(this.dialog.domNode);
		//this.dialog.show();
	}

});

