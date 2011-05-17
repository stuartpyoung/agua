dojo.provide( "plugins.admin.AppRow");

dojo.declare( "plugins.admin.AppRow",
	[ dijit._Widget, dijit._Templated ],
{
	//Path to the template of this widget. 
	templatePath: dojo.moduleUrl("plugins", "admin/templates/approw.html"),

	// Calls dijit._Templated.widgetsInTemplate
	widgetsInTemplate : true,

	// PARENT plugins.admin.Apps WIDGET
	parentWidget : null,

	constructor : function(args)
	{

		this.submitOn = args.submit;
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

		// CONNECT TOGGLE EVENT
		var appRowObject = this;
		dojo.connect( this.name, "onclick", function(event) {
			appRowObject.toggle();
		});

		// ADD 'EDIT' ONCLICKS
		var appRowObject = this;
		var array = [ "executor", "version", "location", "description", "notes" ];
		for ( var i in array )
		{
			dojo.connect(this[array[i]], "onclick", function(event)
				{
					appRowObject.parentWidget.editAppRow(appRowObject, event.target);
					event.stopPropagation(); //Stop Event Bubbling
				}
			);
		}

		// USE this.submitOn TO DECIDE IF CHECKBOX IS SELECTED
		if ( this.submitOn != null && this.submitOn == 1 )
		{

			this.submit.setValue("on");
		}
	},


	// DEBUG	
	submitChange : function ()
	{

		// GET INPUTS
		var inputs = this.parentWidget.getEditedInputs(this);
		if ( inputs == null ) return;

		// SAVE APPLICATION
		this.parentWidget.saveApp(inputs);
	},

	// TOGGLE HIDDEN DETAILS	
	toggle : function ()
	{

		var array = [ "executor", "version", "location", "submitContainer", "description", "notes" ];
		for ( var i in array )
		{

			if ( this[array[i]].style.display == 'table-cell' )	
				this[array[i]].style.display='none';
			else
				this[array[i]].style.display = 'table-cell';
		}
	}

});

