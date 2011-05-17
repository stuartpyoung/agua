dojo.provide( "plugins.admin.AccessRow");


dojo.declare( "plugins.admin.AccessRow",
	[ dijit._Widget, dijit._Templated ],
{
	////////}
//Path to the template of this widget. 
templatePath: dojo.moduleUrl("plugins", "admin/templates/accessrow.html"),

// Calls dijit._Templated.widgetsInTemplate
widgetsInTemplate : true,

// PARENT plugins.admin.Apps WIDGET
parentWidget : null,

// INSTANTIATION ARGUMENTS
args : null,

constructor : function(args) {
	this.args = args;
},

postCreate : function() {

	this.startup();
},

startup : function () {

	this.inherited(arguments);

	for ( var right in this.args )
	{
		if ( right == "groupname" || right == "owner" ) continue;

		if ( this.args[right] == 1 )
		{
			dojo.addClass(this[right], 'allowed');
		}
		else
		{
			dojo.addClass(this[right], 'denied');
		}
	}
},

// TOGGLE HIDDEN NODES
toggle : function (event) {
	var value = dojo.attr(event.target, 'value');
	event.stopPropagation();	

	var ownership = "group";
	if ( value.match(/^world/) )	ownership = "world";

	var permission;
	if ( value.match(/^group/) ) permission = value.replace(/^group/, '');
	if ( value.match(/^world/) ) permission = value.replace(/^world/, '');

	var setClass = "allowed";
	var unsetClass = "denied";
	if ( dojo.hasClass(this[value], "allowed") )
	{
		setClass = "denied";
		unsetClass = "allowed";
	}

	this.descendingOwnership(ownership, permission, setClass, unsetClass);

	if ( ownership == "world" && setClass == "allowed" )
	{
		var deniedOk = false;
		this.descendingOwnership("group", permission, setClass, unsetClass, deniedOk);
	}
},

descendingOwnership : function (ownership, permission, setClass, unsetClass, deniedOk) {

	if ( deniedOk == null )	deniedOk = true;


	if ( permission == "write" )
	{
		this.setClasses(ownership, "write", setClass, unsetClass);
		this.setClasses(ownership, "copy", setClass, unsetClass);	
		this.setClasses(ownership, "view", setClass, unsetClass);
	}

	if ( permission == "copy" )
	{
		this.setClasses(ownership, "view", setClass, unsetClass);
		this.setClasses(ownership, "copy", setClass, unsetClass);
		if ( setClass == "denied" )	
			this.setClasses(ownership, "write", setClass, unsetClass);
		else if ( deniedOk == true )
			this.setClasses(ownership, "write", unsetClass, setClass);
	}

	if ( permission == "view" )
	{
		this.setClasses(ownership, "view", setClass, unsetClass);
		if ( setClass == "denied" && deniedOk == true )	
		{
			this.setClasses(ownership, "write", setClass, unsetClass);	
			this.setClasses(ownership, "copy", setClass, unsetClass);	
		}
		else if ( deniedOk == true )
		{
			this.setClasses(ownership, "write", unsetClass, setClass);	
			this.setClasses(ownership, "copy", unsetClass, setClass);	
		}
	}
},

setClasses : function (ownership, permission, setClass, unsetClass) {

	dojo.addClass(this[ownership + permission ], setClass);
	dojo.removeClass(this[ownership + permission], unsetClass);
}



});

