dojo.provide("plugins.admin.Access");

// ALLOW THE ADMIN USER TO ADD, REMOVE AND MODIFY USERS

// NEW USERS MUST HAVE username AND email

//dojo.require("dijit.dijit"); // optimize: load dijit layer
dojo.require("dijit.form.Button");
dojo.require("dijit.form.TextBox");
dojo.require("dijit.form.Textarea");
dojo.require("dojo.parser");
dojo.require("dojo.dnd.Source");
dojo.require("plugins.core.Common");

// FORM VALIDATION
dojo.require("plugins.form.ValidationTextarea");

// HAS A
dojo.require("plugins.admin.AccessRow");

dojo.declare(
    "plugins.admin.Access",
	[ dijit._Widget, dijit._Templated, plugins.core.Common ],
{

	/////}

//Path to the template of this widget. 
templatePath: dojo.moduleUrl("plugins", "admin/templates/access.html"),

// Calls dijit._Templated.widgetsInTemplate
widgetsInTemplate : true,

// OR USE @import IN HTML TEMPLATE
cssFiles : [ "plugins/admin/css/access.css" ],

// PARENT WIDGET
parentWidget : null,

// ROW WIDGETS
accessRows : new Array,

rights: [ 'groupname', 'groupwrite', 'groupcopy', 'groupview', 'worldwrite', 'worldcopy', 'worldview' ],

constructor : function (args) {


	// GET INFO FROM ARGS
	this.parentWidget = args.parentWidget;
	this.tabContainer = args.tabContainer;

	this.loadCSS();
},

postMixInProperties: function() {
},

postCreate: function() {
	this.startup();
},

startup : function () {

	// COMPLETE CONSTRUCTION OF OBJECT
	this.inherited(arguments);	 

	// ADD ADMIN TAB TO TAB CONTAINER		
	this.tabContainer.addChild(this.accessTab);
	this.tabContainer.selectChild(this.accessTab);

	this.buildTable();
},

// BUILD USERS TABLE
buildTable : function () {


	// GET ACCESS TABLE DATA
	var accessArray = Agua.getAccess();

	// CLEAN TABLE
	if ( this.table.childNodes )
	{
		while ( this.table.childNodes.length > 2 )
		{
			this.table.removeChild(this.table.childNodes[2]);
		}
	}

	// BUILD ROWS
	this.tableRows = [];
	for ( var rowCounter = 0; rowCounter < accessArray.length; rowCounter++)
	{
		var accessRow = new plugins.admin.AccessRow(accessArray[rowCounter]);
		this.table.appendChild(accessRow.row);
		this.tableRows.push(accessRow);
	}

},	// buildTable

togglePermission : function (event) {
	var parentNode = event.target.parentNode;

	var nodeClass = event.target.getAttribute('class');
	if ( nodeClass.match('allowed') )
	{
		dojo.removeClass(event.target,'allowed');
		dojo.addClass(event.target,'denied');
	}
	else
	{
		dojo.removeClass(event.target,'denied');
		dojo.addClass(event.target,'allowed');
	}
},

saveAccess : function () {

	// COLLECT DATA HERE
	var dataArray = new Array;
	for ( var i = 0; i < this.tableRows.length; i++ )
	{
		var data = new Object;
		data.owner 		= Agua.cookie('username');
		data.groupname	= this.tableRows[i].args.groupname;
		for ( var j = 1; j < this.rights.length; j++ )
		{
			data[this.rights[j]] = 0;
			if ( dojo.hasClass(this.tableRows[i][this.rights[j]], 'allowed') )
				data[this.rights[j]] = 1;
		}
		dataArray.push(data);
	}

	// CREATE JSON QUERY
	var url = Agua.cgiUrl + "agua?";
	var query = new Object;
	query.username = Agua.cookie('username');
	query.sessionId = Agua.cookie('sessionId');
	query.mode = "saveAccess";
	query.data = dataArray;


	// SEND TO SERVER
	dojo.xhrPut(
		{
			url: url,
			contentType: "text",
			putData: dojo.toJson(query),
			timeout: 3000,
			load: function(response, ioArgs) {
				return response;
			},
			error: function(response, ioArgs) {
				return response;
			}
		}
	);
}


}); // plugins.admin.Access

