dojo.provide("plugins.workflow.AppsMenu");

// ALLOW THE USER TO ADD, REMOVE AND MODIFY APPS

//dojo.require("dijit.dijit"); // optimize: load dijit layer
dojo.require("dojo.parser");

// HAS A
dojo.require("dijit.Menu");

// INHERITS
dojo.require("plugins.core.Common");


dojo.declare(
    "plugins.workflow.AppsMenu",
	[ dijit._Widget, dijit._Templated, plugins.core.Common ],
{

//Path to the template of this widget. 
templatePath: dojo.moduleUrl("plugins", "workflow/templates/appsmenu.html"),

// Calls dijit._Templated.widgetsInTemplate
widgetsInTemplate : true,

//addingApp STATE
addingApp : false,

// OR USE @import IN HTML TEMPLATE
cssFiles : [ dojo.moduleUrl("plugins") + "/workflow/css/appsmenu.css" ],


// PARENT WIDGET
parentWidget : null,

// CORE WORKFLOW OBJECTS
core : null,

/////}
constructor : function(args)
{
	// GET INFO FROM ARGS
	this.parentWidget = args.parentWidget;

	// LOAD SORIA AND FILEPICKER CSS
	this.loadCSS();		
},

postCreate : function()
{

	this.startup();
},


startup : function ()
{

	// COMPLETE CONSTRUCTION OF OBJECT
	this.inherited(arguments);	 

	// SET DRAG APP - LIST OF APPS
	this.setMenu();
},


// BIND THE MENU TO A NODE
bind : function (node)
{

	if ( node == null )
	{

	}
	return this.menu.bindDomNode(node);	
},

// SHOW 'ABOUT' INFORMATION
about : function (event)
{

		event.stopPropagation();

},		

// OPEN WINDOW TO APPLICATION WEBSITE
website : function (event)
{
	event.stopPropagation();


},

// ADD PROGRAMMATIC CONTEXT MENU
setMenu : function ()
{

	dojo.connect(this.aboutNode, "onClick", dojo.hitch(this, function(event)
	{
		this.about(event);
		event.stopPropagation();
	}));

	dojo.connect(this.websiteNode, "onClick", dojo.hitch(this, function(event)
	{
		this.website(event);
		event.stopPropagation();
	}));	
}

}); // plugins.workflow.AppsMenu

