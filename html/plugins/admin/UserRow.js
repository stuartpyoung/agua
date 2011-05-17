dojo.provide( "plugins.admin.UserRow");


dojo.declare( "plugins.admin.UserRow",
	[ dijit._Widget, dijit._Templated ],
{
	//Path to the template of this widget. 
	templatePath: dojo.moduleUrl("plugins", "admin/templates/userrow.html"),

	// Calls dijit._Templated.widgetsInTemplate
	widgetsInTemplate : true,

	// PARENT plugins.admin.Users WIDGET
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

		//dojo.connect( this.username, "onclick", this.toggle);
		var userRowObject = this;
		dojo.connect( this.username, "onclick", function(event) {

			userRowObject.toggle();

			event.stopPropagation(); //Stop Event Bubbling 			
		});

		// ADD 'EDIT' ONCLICK
		var userRowObject = this;
		dojo.connect(this.firstname, "onclick", function(event)
			{

				userRowObject.parentWidget.editUserRow(userRowObject, event.target);
				event.stopPropagation(); //Stop Event Bubbling 			
			}
		);

		// ADD 'EDIT' ONCLICK
		var userRowObject = this;
		dojo.connect(this.lastname, "onclick", function(event)
			{

				userRowObject.parentWidget.editUserRow(userRowObject, event.target);
				event.stopPropagation(); //Stop Event Bubbling 			
			}
		);

		// ADD 'EDIT' ONCLICK
		var userRowObject = this;
		dojo.connect(this.email, "onclick", function(event)
			{

				userRowObject.parentWidget.editUserRow(userRowObject, event.target);
				event.stopPropagation(); //Stop Event Bubbling 			
			}
		);

		// ADD 'PASSWORD' ONCLICK
		var userRowObject = this;
		dojo.connect(this.password, "onclick", function(event)
			{

				userRowObject.clearValue(event.target);
				userRowObject.parentWidget.editUserRow(userRowObject, event.target);
				event.stopPropagation(); //Stop Event Bubbling 			
			}
		);

	},

	toggle : function ()
	{

		if ( this.email.style.display == 'block' ) this.email.style.display='none';
		else this.email.style.display = 'block';
		if ( this.firstname.style.display == 'inline' ) this.firstname.style.display='none';
		else this.firstname.style.display = 'inline';
		if ( this.lastname.style.display == 'inline' ) this.lastname.style.display='none';
		else this.lastname.style.display = 'inline';
		if ( this.password.style.display == 'block' ) this.password.style.display='none';
		else this.password.style.display = 'block';
	},

	clearValue : function (element) {

		if ( element.clicked == true ) return;

		element.clicked = true;
		element.innerHTML = '';
		element.focus();
	}
});
