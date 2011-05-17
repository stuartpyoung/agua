dojo.provide("plugins.workflow.HistoryPane");

// DISPLAY THE STATUS OF A WORKFLOW STAGE

dojo.declare( "plugins.workflow.HistoryPane",
	[ dijit._Widget, dijit._Templated ],
{
//Path to the template of this widget. 
templatePath: dojo.moduleUrl("plugins", "workflow/templates/historypane.html"),

// Calls dijit._Templated.widgetsInTemplate
widgetsInTemplate : true,

srcNodeRef: null,

// ROWS OF ENTRIES
rows : null,

// CORE WORKFLOW OBJECTS
core : null,

/////}
constructor : function(args)
{		
	//this.project = rows[0].project;
	//this.workflow = rows[0].workflow;
	this.rows = args.rows;

	this.core = args.core;

},

postMixInProperties: function()
{
},

postCreate: function()
{

	this.startup();		
},

// startup
//
startup : function ()
{

	// SET UP THE ELEMENT OBJECTS AND THEIR VALUE FUNCTIONS
	this.inherited(arguments);

	for ( var i = 0; i < this.rows.length; i++ )
	{
		var historyPaneRow = new plugins.workflow.HistoryPaneRow(this.rows[i]);

		this.rowsNode.innerHTML += historyPaneRow.domNode.innerHTML;

	}
}



}); // end of plugins.workflow.HistoryPane



dojo.declare( "plugins.workflow.HistoryPaneRow",
	[ dijit._Widget, dijit._Templated ],
{

//Path to the template of this widget. 
templatePath: dojo.moduleUrl("plugins", "workflow/templates/historypanerow.html"),

// Calls dijit._Templated.widgetsInTemplate
widgetsInTemplate : true,

srcNodeRef: null,

// CORE WORKFLOW OBJECTS
core : null,

/////}
constructor : function(args)
{		
	this.core = args.core;
},

// postMixInProperties
//Inherited from dijit._Widget and called just before template
//instantiation in buildRendering. This method is especially useful
//for manipulating the template before it becomes visible.
postMixInProperties: function()
{
},

postCreate : function()
{

	this.startup();
},

startup : function ()
{

	this.inherited(arguments);


}

});
