dojo.provide( "plugins.workflow.AppType");


dojo.declare( "plugins.workflow.AppType",
	[ dijit._Widget, dijit._Templated ],
{
	//Path to the template of this widget. 
	templatePath: dojo.moduleUrl("plugins", "workflow/templates/apptype.html"),

	// Calls dijit._Templated.widgetsInTemplate
	widgetsInTemplate : true,

	// PARENT plugins.workflow.Apps WIDGET
	parentWidget : null,

	constructor : function(args)
	{
		this.parentWidget = args.parentWidget;
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
