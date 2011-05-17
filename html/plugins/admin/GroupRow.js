dojo.provide( "plugins.admin.GroupRow");


dojo.declare( "plugins.admin.GroupRow",
	[ dijit._Widget, dijit._Templated ],
{
	//Path to the template of this widget. 
	templatePath: dojo.moduleUrl("plugins", "admin/templates/grouprow.html"),

	// Calls dijit._Templated.widgetsInTemplate
	widgetsInTemplate : true,

	// PARENT plugins.admin.Groups WIDGET
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

		//dojo.connect( this.name, "onclick", this.toggle);
		var groupRowObject = this;
		dojo.connect( this.name, "onclick", function(event) {

			groupRowObject.toggle();

			event.stopPropagation(); //Stop Event Bubbling 			
		});

		// ADD 'EDIT' ONCLICK
		var groupRowObject = this;

		// DESCRIPTION
		dojo.connect(this.description, "onclick", function(event)
			{
				groupRowObject.parentWidget.editGroupRow(groupRowObject, event.target);
				//event.stopPropagation(); //Stop Event Bubbling
			}
		);

		// NOTES
		dojo.connect(this.notes, "onclick", function(event)
			{
				groupRowObject.parentWidget.editGroupRow(groupRowObject, event.target);
				//event.stopPropagation(); //Stop Event Bubbling 			
			}
		);
	},

	toggle : function ()
	{

		if ( this.description.style.display == 'block' ) this.description.style.display='none';
		else this.description.style.display = 'block';
		if ( this.notes.style.display == 'block' ) this.notes.style.display='none';
		else this.notes.style.display = 'block';
	}

});

