dojo.provide( "plugins.admin.GroupProjectRow");


dojo.declare( "plugins.admin.GroupProjectRow",
	[ dijit._Widget, dijit._Templated ],
{
	//Path to the template of this widget. 
	templatePath: dojo.moduleUrl("plugins", "admin/templates/groupprojectrow.html"),

	// Calls dijit._Templated.widgetsInTemplate
	widgetsInTemplate : true,

	// PARENT plugins.admin.Sources WIDGET
	parentWidget : null,

	constructor : function(args)
	{

		this.parentWidget = args.parentWidget;
		//this.inherited(arguments);
	},

	postCreate : function()
	{

		this.startup();
	},

	startup : function ()
	{

		this.inherited(arguments);

		var groupProjectRowObject = this;
		dojo.connect( this.name, "onclick", function(event) {

			groupProjectRowObject.toggle();
			event.stopPropagation(); //Stop Event Bubbling 			
		});

		//// ADD 'EDIT' ONCLICK
		//var groupProjectRowObject = this;
		//dojo.connect(this.description, "onclick", function(event)
		//	{
		//
		//		groupProjectRowObject.parentWidget.editGroupProjectRow(groupProjectRowObject, event.target);
		//		event.stopPropagation(); //Stop Event Bubbling 			
		//	}
		//);
		//
		//// ADD 'EDIT' ONCLICK
		//var groupProjectRowObject = this;
		//dojo.connect(this.location, "onclick", function(event)
		//	{
		//
		//		groupProjectRowObject.parentWidget.editGroupProjectRow(groupProjectRowObject, event.target);
		//		event.stopPropagation(); //Stop Event Bubbling 			
		//	}
		//);
	},

	toggle : function ()
	{

		if ( this.description.style.display == 'block' ) this.description.style.display='none';
		else this.description.style.display = 'block';
	}
});
	