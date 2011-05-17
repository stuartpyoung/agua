dojo.provide( "plugins.admin.SourceRow");


dojo.declare( "plugins.admin.SourceRow",
	[ dijit._Widget, dijit._Templated ],
{
	//Path to the template of this widget. 
	templatePath: dojo.moduleUrl("plugins", "admin/templates/sourcerow.html"),

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

		//dojo.connect( this.name, "onclick", this.toggle);
		var sourceRowObject = this;
		dojo.connect( this.name, "onclick", function(event) {

			sourceRowObject.toggle();

			event.stopPropagation(); //Stop Event Bubbling 			
		});

		// ADD 'EDIT' ONCLICK
		var sourceRowObject = this;
		dojo.connect(this.description, "onclick", function(event)
			{

				sourceRowObject.parentWidget.editSourceRow(sourceRowObject, event.target);
				event.stopPropagation(); //Stop Event Bubbling 			
			}
		);

		// ADD 'EDIT' ONCLICK
		var sourceRowObject = this;
		dojo.connect(this.location, "onclick", function(event)
			{

				sourceRowObject.parentWidget.editSourceRow(sourceRowObject, event.target);
				event.stopPropagation(); //Stop Event Bubbling 			
			}
		);
	},

	toggle : function ()
	{

		if ( this.description.style.display == 'block' ) this.description.style.display='none';
		else this.description.style.display = 'block';
		if ( this.location.style.display == 'block' ) this.location.style.display='none';
		else this.location.style.display = 'block';
	}

});
