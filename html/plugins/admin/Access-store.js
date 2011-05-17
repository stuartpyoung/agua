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
dojo.require("plugins.admin.UserRow");

dojo.declare(
    "plugins.admin.Access",
	[ dijit._Widget, dijit._Templated, plugins.core.Common ],
{
	//Path to the template of this widget. 
	templatePath: dojo.moduleUrl("plugins", "admin/templates/access.html"),

	// Calls dijit._Templated.widgetsInTemplate
	widgetsInTemplate : true,

	// Calls dijit._Templated.widgetsInTemplate
	widgetsInTemplate : true,

	//title : "<span class='adminHeader'> Access </span> Manage Access to Groups",        
	objectStore : null,
	tableId : "adminAccessTable",
	tableRowConnections: new Array,
	deleteRadioId : "adminAccessDeleteRadio",

	// OR USE @import IN HTML TEMPLATE
	cssFiles : [ "plugins/admin/css/access.css" ],

	// PARENT WIDGET
	parentWidget : null,

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

	groupTable : function () {

	},
	// BUILD USERS TABLE
	buildTable : function () {

		// GET ACCESS TABLE DATA
		var accessArray = Agua.getAccess();

		// SET STORE
		var data = {identifier: "project", items: []};
		for ( var i = 0; i < accessArray.length; i++ )
		{
			data.items[i] = accessArray[i];
		}

		var objectStore = new dojo.data.ItemFileWriteStore(	{	data: data	}	);
		this.objectStore = objectStore;

		//General Structure
		// http://dojotoolkit.org/reference-guide/dojo/data/ItemFileReadStore.html		
		//The ItemFileReadStore expects a specific structure to its data, as defined below:		
		//{
		//	"label": "some attribute",   //Optional attribute used to indicate which attribute on an item should act as a human-readable label for display purposes.
		//	
		//	
		//	"identifier": "some attribute",  //Optional attribute used to indicate which attribute on an item acts as a unique identifier for that item. If it is not defined, then the ItemFileReadStore will simply number the items and use that number as a unique index to the item.
		//	
		//	
		//	"items:" [  //The array of JavaScript objects that act as the root items of the data store
		//			{ /* Some set of name/value attributes */ },
		//			{ /* ... */ },
		//			...
		//	]
		//}

		// WRITE TABLE BASED ON ITEMS IN STORE REQUEST
		var adminAccess = this;
		var groupTable = function(items, request)
		{

			// CREATE TABLE adminAccessTable 
			if ( document.getElementById(adminAccess.tableId) )
			{
				adminAccess.domNode.removeChild(document.getElementById(adminAccess.tableId));
			}

			var table = document.createElement("table");
			table.setAttribute('class', 'adminTable');
			table.id = adminAccess.tableId;


			// APPEND TABLE TO RIGHT PANE

			adminAccess.domNode.appendChild(table);

			// CREATE HEADER ROW
			var headerRow = document.createElement("tr");
			headerRow.setAttribute('class', 'adminHeaderRow');
			var headers = [ "Group", "Project", "Owner rights", "Group rights" , "World rights" ];
			var cols = [ 1, 1, 3, 3, 3 ];
			for ( var i = 0; i < headers.length; i++ )
			{
				var tableData = document.createElement("td");
				tableData.setAttribute('colspan', cols[i]);
				//tableData.setAttribute('class', 'adminAccessTableHeader');
				var text = document.createTextNode(headers[i]);
				tableData.appendChild(text);
				headerRow.appendChild(tableData);
			}
			table.appendChild(headerRow);	

			var subHeaderRow = document.createElement("tr");
			subHeaderRow.setAttribute('class', 'adminHeaderRow');
			var subHeaders = [ "", "", "Edit", "Copy", "View", "Edit", "Copy", "View", "Edit", "Copy", "View" ];
			for ( var i = 0; i < subHeaders.length; i++ )
			{
				var tableData = document.createElement("td");
				var text = document.createTextNode(subHeaders[i]);
				tableData.setAttribute('class', 'adminAccessTableSubHeader');
				tableData.appendChild(text);
				subHeaderRow.appendChild(tableData);
			}
			table.appendChild(subHeaderRow);	


			// SHOW PROJECT ACCESS INFO
			for ( var rowCounter = 0; rowCounter < items.length; rowCounter++)
			{
				// CREATE TABLE ROW FOR EACH USER
				var tableRow = document.createElement("tr");
				tableRow.setAttribute('class', 'adminAccessTableRow');
				tableRow.setAttribute('id', 'adminAccessTableRow' + rowCounter);
				var item = items[rowCounter];
				var data = [ adminAccess.objectStore.getValue(item, 'groupname'), adminAccess.objectStore.getLabel(item), adminAccess.objectStore.getValue(item, 'ownerrights'), adminAccess.objectStore.getValue(item, 'grouprights'), adminAccess.objectStore.getValue(item, 'worldrights')];

				var columns = [];
				columns.push(data[0]);
				columns.push(data[1]);
				for ( var i = 2; i < data.length; i++ )
				{
					if ( data[i] == 7 )	{	columns.push(true, true, true);		}
					if ( data[i] == 6 )	{	columns.push(true, true, false);	}
					if ( data[i] == 5 )	{	columns.push(true, false, true);	}
					if ( data[i] == 4 )	{	columns.push(true, false, false);	}
					if ( data[i] == 3 )	{	columns.push(false, true, true);	}
					if ( data[i] == 2 )	{	columns.push(false, true, false);	}
					if ( data[i] == 1 )	{	columns.push(false, false, true);	}
					if ( data[i] == 0 )	{	columns.push(false, false, false);	}
				}


				// SHOW COLUMNS: USERNAME, FULL NAME, PASSWORD
				for ( var columnCounter = 0; columnCounter < columns.length; columnCounter++ )
				{
					var tableData = document.createElement("td");

					if ( columnCounter == 0 || columnCounter == 1 )
					{
						tableData.setAttribute('class', 'adminAccessTableData');	
						var text = document.createTextNode(columns[columnCounter]);
						tableData.appendChild(text);
					}
					else
					{
						if ( columns[columnCounter] == true )
						{
							tableData.setAttribute('class', 'adminAccessAllowed');	
						}
						else
						{
							tableData.setAttribute('class', 'adminAccessDenied');	
						}

						// ADD ONCLICK TO TOGGLE CLASS
						dojo.connect(tableData, "onclick", function(event){
							adminAccess.togglePermission(event);
						});
					}

					// APPEND TABLE DATA TO TABLE ROW
					tableRow.appendChild(tableData);							
				}

				// APPEND TABLE DATA TO TABLE ROW
				tableRow.appendChild(tableData);			

				// APPEND ROW TO TABLE
				table.appendChild(tableRow);
			}
		}; // groupTable

		var returned = this.objectStore.fetch( { query: "*", onComplete: groupTable});
		//var returned = this.objectStore.fetch( { query: {}, onComplete: groupTable});

	},	// buildTable


	togglePermission : function (event) {

		var target = event.target;

		var nodeClass = target.getAttribute('class');
		if ( nodeClass == 'adminAccessAllowed' )
		{
			target.setAttribute('class', 'adminAccessDenied');
		}
		else
		{
			target.setAttribute('class', 'adminAccessAllowed');
		}
	},


	//// MAKE SURE THAT PERMISSIONS ARE NUMERIC AND BETWEEN 1 AND 7
	//checkPermissions : function (text)
	//{
	//	if ( ! text.match(/^\d+$/) ) return 1;
	//	if ( text < 1 )	return 1;
	//	if ( text > 7 ) return 1;
	//	
	//	return text;
	//},

	// CURRENTLY UNUSED CALLBACK FOR saveRow
	updateAccess : function (json) {
	},




	saveStore : function () {

		// COLLECT DATA HERE
		var dataArray = new Array;

		var rights = [ "ownerrights", "grouprights", "worldrights" ];
		var table = dojo.byId(this.tableId);
		var rows = table.childNodes;
		for ( var i = 2; i < rows.length; i++ )
		{
			var data = new Object;
			var tableDatas = rows[i].childNodes;
			data.project = tableDatas[0].innerHTML.toString();
			var rightsCounter = 0;
			for ( var j = 1; j < tableDatas.length; j++ )
			{
				var value = 0;
				var values = [ 4, 2, 1];
				for ( var k = 0; k < 3; k++ )
				{
					var node = tableDatas[(j + k)];
					var nodeClass = node.getAttribute('class').toString();
					if ( nodeClass == 'adminAccessAllowed' )
					{
						value += values[k];
					}
				}

				data[rights[rightsCounter]] = value;
				rightsCounter++;

				j+=2;
			}

			dataArray.push(data);
		}



		// SAVE FOR LATER: DO IT USING THE STORE			
		////////// print out for debugging
		////////var adminAccess = this;
		////////
		////////var store = this.objectStore;
		////////
		////////for ( var i = 0; i < store._arrayOfTopLevelItems.length; i++ )
		////////{
		////////}
		////////
		////////// SEND TO SERVER AS RAW XHR POST
		////////var dataArray = new Array;
		////////for ( var i = 0; i < store._arrayOfTopLevelItems.length; i++ )
		////////{
		////////	
		////////	var data = new Object;
		////////	var item = store._arrayOfTopLevelItems[i];
		////////	var attributes = store.getAttributes(item);
		////////	if (attributes && attributes.length > 0)
		////////	{
		////////		for ( var j = 0; j < attributes.length; j++ )
		////////		//for ( var j = 0; j < 1; j++ )
		////////		{
		////////			var values = store.getValues(item, attributes[j]);
		////////
		////////			if ( values )
		////////			{
		////////				
		////////				// MULTI-VALUE ATTRIBUTE
		////////				if (values.length > 1 )
		////////				{
		////////					data[attributes[i]] = [];
		////////					for ( var k = 0; k < values.length; k++ )
		////////					{
		////////						var value = values[k];
		////////						data[attributes[j]].push(value);
		////////					}
		////////				}
		////////				// SINGLE VALUE ATTRIBUTE
		////////				else
		////////				{
		////////					data[attributes[j]] = values[0];
		////////				}
		////////			}
		////////		}
		////////		
		////////	}
		////////	dataArray.push(data);
		////////}

		var url = Agua.cgiUrl + "agua?";

		// CREATE JSON QUERY
		var query = new Object;
		query.username = Agua.cookie('username');
		query.sessionId = Agua.cookie('sessionId');
		query.mode = "saveAccess";
		query.data = dojo.toJson(dataArray);

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

