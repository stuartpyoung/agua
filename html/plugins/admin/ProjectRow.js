dojo.provide( "plugins.admin.ProjectRow");


dojo.declare( "plugins.admin.ProjectRow",
	[ dijit._Widget, dijit._Templated ],
{
	//Path to the template of this widget. 
	templatePath: dojo.moduleUrl("plugins", "admin/templates/projectrow.html"),

	// Calls dijit._Templated.widgetsInTemplate
	widgetsInTemplate : true,

	// PARENT plugins.admin.Projects WIDGET
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
		var projectRowObject = this;
		dojo.connect( this.name, "onclick", function(event) {

			projectRowObject.toggle();

			event.stopPropagation(); //Stop Event Bubbling 			
		});

		// ADD 'EDIT' ONCLICK
		var projectRowObject = this;
		dojo.connect(this.description, "onclick", function(event)
			{

				projectRowObject.parentWidget.editProjectRow(projectRowObject, event.target);
				event.stopPropagation(); //Stop Event Bubbling
			}
		);

		// ADD 'EDIT' ONCLICK
		var projectRowObject = this;
		dojo.connect(this.notes, "onclick", function(event)
			{

				projectRowObject.parentWidget.editProjectRow(projectRowObject, event.target);
				event.stopPropagation(); //Stop Event Bubbling 			
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

