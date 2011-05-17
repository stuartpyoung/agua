dojo.provide( "plugins.admin.GroupUserRow");


dojo.declare( "plugins.admin.GroupUserRow",
	[ dijit._Widget, dijit._Templated ],
{
	//Path to the template of this widget. 
	templatePath: dojo.moduleUrl("plugins", "admin/templates/groupuserrow.html"),

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

		var groupUserRowObject = this;
		dojo.connect( this.username, "onclick", function(event) {
			groupUserRowObject.toggle();
			event.stopPropagation(); //Stop Event Bubbling 			
		});

		//// ADD 'EDIT' ONCLICK
		//var groupUserRowObject = this;
		//dojo.connect(this.description, "onclick", function(event)
		//	{
		//
		//		groupUserRowObject.parentWidget.editGroupUserRow(groupUserRowObject, event.target);
		//		event.stopPropagation(); //Stop Event Bubbling 			
		//	}
		//);
		//
		//// ADD 'EDIT' ONCLICK
		//var groupUserRowObject = this;
		//dojo.connect(this.location, "onclick", function(event)
		//	{
		//
		//		groupUserRowObject.parentWidget.editGroupUserRow(groupUserRowObject, event.target);
		//		event.stopPropagation(); //Stop Event Bubbling 			
		//	}
		//);
	},

	toggle : function ()
	{

		if ( this.email.style.display == 'block' ) this.email.style.display='none';
		else this.email.style.display = 'block';
		if ( this.location.style.display == 'block' ) this.location.style.display='none';
		else this.location.style.display = 'block';
		if ( this.fullname.style.display == 'block' ) this.fullname.style.display='none';
		else this.fullname.style.display = 'block';
		if ( this.description.style.display == 'block' ) this.description.style.display='none';
		else this.description.style.display = 'block';
	}
});

