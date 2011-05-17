dojo.provide("plugins.workflow.Status");

// DISPLAY THE STATUS OF A WORKFLOW STAGE

dojo.declare( "plugins.workflow.Status",
	[ dijit._Widget, dijit._Templated, plugins.core.Common ],
	//[ plugins.core.Common, plugins.core.Template ],
{
//Path to the template of this widget. 
templatePath: dojo.moduleUrl("plugins", "workflow/templates/status.html"),

// Calls dijit._Templated.widgetsInTemplate
widgetsInTemplate : true,

srcNodeRef: null,

// CORE WORKFLOW OBJECTS
core : null,

/////}

// CONSTRUCTOR	

// Any initialization code would go here in the constructor.
// plugins.report.Template and its superclasses dijit._Widget and
// dijit._Templated do not have parameters in their constructors, so
// there wouldn't be any multiple-inheritance complications
// if you were to include some paramters here.
constructor : function(args)
{		
	this.core = args.core;

	// SET ARGS
	this.name = args.name;
	this.number = args.number;
	this.status = args.status;
	this.started = args.started;
	this.completed = args.completed;
},

//Inherited from dijit._Widget and called just before template
//instantiation in buildRendering. This method is especially useful
//for manipulating the template before it becomes visible.
postMixInProperties: function()
{
},

//You can override this method to manipulate widget once it is
//placed in the UI, but be warned that any child widgets contained
//in it may not be ready yet.        
postCreate: function()
{

},


// startup
//
//
startup : function ()
{

	// SET UP THE ELEMENT OBJECTS AND THEIR VALUE FUNCTIONS
	this.inherited(arguments);

}



}); // end of plugins.workflow.Status
