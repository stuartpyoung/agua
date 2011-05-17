dojo.provide( "plugins.admin.GroupSourceRow");


dojo.declare( "plugins.admin.GroupSourceRow",
	[ dijit._Widget, dijit._Templated ],
{
	//Path to the template of this widget. 
	templatePath: dojo.moduleUrl("plugins", "admin/templates/groupsourcerow.html"),

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

		var groupSourceRowObject = this;
		dojo.connect( this.name, "onclick", function(event) {
			groupSourceRowObject.toggle();
			event.stopPropagation(); //Stop Event Bubbling 			
		});

		//// ADD 'EDIT' ONCLICK
		//var groupSourceRowObject = this;
		//dojo.connect(this.description, "onclick", function(event)
		//	{
		//
		//		groupSourceRowObject.parentWidget.editGroupSourceRow(groupSourceRowObject, event.target);
		//		event.stopPropagation(); //Stop Event Bubbling 			
		//	}
		//);
		//
		//// ADD 'EDIT' ONCLICK
		//var groupSourceRowObject = this;
		//dojo.connect(this.location, "onclick", function(event)
		//	{
		//
		//		groupSourceRowObject.parentWidget.editGroupSourceRow(groupSourceRowObject, event.target);
		//		event.stopPropagation(); //Stop Event Bubbling 			
		//	}
		//);
	},

	toggle : function ()
	{

		//if ( this.location.style.display == 'block' ) this.location.style.display='none';
		//else this.location.style.display = 'block';
		if ( this.description.style.display == 'block' ) this.description.style.display='none';
		else this.description.style.display = 'block';
	}
});

