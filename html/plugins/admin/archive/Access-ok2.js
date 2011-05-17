dojo.provide("plugins.admin.Access");

dojo.require("dijit.dijit"); // optimize: load dijit layer	[ dijit._Widget, dijit._Templated, plugins.core.Common ],

//dojo.require("dijit.dijit"); // optimize: load dijit layer
dojo.require("dijit.form.Textarea");
dojo.require("dojo.dnd.Source");


dojo.require("dijit.form.Button");
dojo.require("dijit.form.TextBox");
dojo.require("dijit.form.ValidationTextBox");
dojo.require("dijit.form.NumberTextBox");
dojo.require("dojo.data.ItemFileWriteStore");
dojo.require("dojo.parser");
dojo.require("plugins.core.Common");


//http://www.ensembl.org/biomart/martview/e61a39e4e8e306354f8f2a7a70bbc53c/e61a39e4e8e306354f8f2a7a70bbc53c
//http://localhost:8080/Bioptic0.2.5/html/dojo-1.5.0/demos/
//http://localhost:8080/Bioptic0.2.5/html/dojo-1.5.0/dojox/form/tests/test_SelectStack.html
//Real World Dojo part 5: Custom Components
//http://www.dotnetmafia.com/blogs/jamesashley/archive/2008/10/28/761.aspx


dojo.declare(
    "plugins.admin.Access",
	[ dijit._Widget, dijit._Templated, plugins.core.Common ],
{
	//Path to the template of this widget. 
	templatePath: dojo.moduleUrl("plugins", "admin/templates/access.html"),

	//// Calls dijit._Templated.widgetsInTemplate
	//widgetsInTemplate : true,
	//
	////title : "<span class='adminHeader'> Access </span> Manage Access to Groups",        
	//url : '',
	//id : '',
	//filename: '',
	//loading : null,
	//rightPaneId : '',
	//groupStore : null,
	//groupCombo : null,
	//tableId : "adminAccessTable",
	//tableRowConnections: new Array,
	//deleteRadioId : "adminAccessDeleteRadio",
	//
	//// OR USE @import IN HTML TEMPLATE
	//cssFiles : [ "plugins/admin/css/access.css" ],
	//
	//// PARENT WIDGET
	//parentWidget : null,

//	// Any initialization code would go here in the constructor.
//	// plugins.report.Template and its superclasses dijit._Widget and
//	// dijit._Templated do not have parameters in their constructors, so
//	// there wouldn't be any multiple-inheritance complications
//	// if you were to include some paramters here.
//	constructor : function (args)
//	{
//		//plugins.admin.Access.superclass.constructor (args);
//		
//
//		// GET INFO FROM ARGS
//		this.parentWidget = args.parentWidget;
//		this.tabContainer = args.tabContainer;
//
//		//this.rightPaneId = args.rightPaneId;
//
//		this.loadCSS();
//	},
//
//	//Inherited from dijit._Widget and called just before template
//	//instantiation in buildRendering. This method is especially useful
//	//for manipulating the template before it becomes visible.
//	postMixInProperties: function()
//	{
//		//this.popup = new dijit.Dialog({});
//	},
//
//
//	//You can override this method to manipulate widget once it is
//	//placed in the UI, but be warned that any child widgets contained
//	//in it may not be ready yet.        
//	postCreate: function()
//	{
//		this.startup();
//	},
//
//
//	startup : function ()
//	{
//
//		// COMPLETE CONSTRUCTION OF OBJECT
//		this.inherited(arguments);	 
//
//		// ADD ADMIN TAB TO TAB CONTAINER		
//		this.tabContainer.addChild(this.accessTab);
//		this.tabContainer.selectChild(this.accessTab);
//
//		this.buildTable();
//
//	},
//
//	// BUILD USERS TABLE
//	buildTable : function ()
//	{
//		
//		// GET GROUP STORE JSON
//		var adminAccess = this;
//		
//		var url = Agua.cgiUrl + "agua?";
//		var query = new Object;
//		query.username = Agua.cookie('username');
//		query.sessionId = Agua.cookie('sessionId');
//		query.mode = "getAccess";
//		if ( ! this.groupStore )
//		{
//			// SEND TO SERVER
//			dojo.xhrPut(
//				{
//					url: url,
//					contentType: "text",
//					sync : true,
//					handleAs: "json",
//					putData: dojo.toJson(query),
//					timeout: 3000,
//					load: function(data)
//					{
//						var groupStore = new dojo.data.ItemFileWriteStore(	{	data: data	}	);
//						adminAccess.groupStore = groupStore;
//					},
//					error: function(response, ioArgs) {
//					}
//				}
//			);
//		}
//
//
//		// WRITE TABLE BASED ON ITEMS IN STORE REQUEST
//		var groupTable = function(items, request)
//		{
//			
//			// CREATE TABLE adminAccessTable 
//			if ( document.getElementById(adminAccess.tableId) )
//			{
//				adminAccess.domNode.removeChild(document.getElementById(adminAccess.tableId));
//			}
//			
//			var table = document.createElement("table");
//			table.setAttribute('class', 'adminTable');
//			table.id = adminAccess.tableId;
//
//				
//			// APPEND TABLE TO RIGHT PANE
//			adminAccess.domNode.appendChild(table);
//
//			// CREATE HEADER ROW
//			var headerRow = document.createElement("tr");
//			headerRow.setAttribute('class', 'adminHeaderRow');
//			var headers = [ "Group", "Project", "Owner rights", "Group rights" , "World rights" ];
//			var cols = [ 1, 1, 3, 3, 3 ];
//			for ( var i = 0; i < headers.length; i++ )
//			{
//				var tableData = document.createElement("td");
//				tableData.setAttribute('colspan', cols[i]);
//				//tableData.setAttribute('class', 'adminAccessTableHeader');
//				var text = document.createTextNode(headers[i]);
//				tableData.appendChild(text);
//				headerRow.appendChild(tableData);
//			}
//			table.appendChild(headerRow);	
//
//			var subHeaderRow = document.createElement("tr");
//			subHeaderRow.setAttribute('class', 'adminHeaderRow');
//			var subHeaders = [ "", "", "Edit", "Copy", "View", "Edit", "Copy", "View", "Edit", "Copy", "View" ];
//			for ( var i = 0; i < subHeaders.length; i++ )
//			{
//				var tableData = document.createElement("td");
//				var text = document.createTextNode(subHeaders[i]);
//				tableData.setAttribute('class', 'adminAccessTableSubHeader');
//				tableData.appendChild(text);
//				subHeaderRow.appendChild(tableData);
//			}
//			table.appendChild(subHeaderRow);	
//
//
//			// SHOW PROJECT ACCESS INFO
//			for ( var rowCounter = 0; rowCounter < items.length; rowCounter++)
//			{
//				// CREATE TABLE ROW FOR EACH USER
//				var tableRow = document.createElement("tr");
//				tableRow.setAttribute('class', 'adminAccessTableRow');
//				tableRow.setAttribute('id', 'adminAccessTableRow' + rowCounter);
//				var item = items[rowCounter];
//				var data = [ adminAccess.groupStore.getValue(item, 'groupname'), adminAccess.groupStore.getLabel(item), adminAccess.groupStore.getValue(item, 'ownerrights'), adminAccess.groupStore.getValue(item, 'grouprights'), adminAccess.groupStore.getValue(item, 'worldrights')];
//				
//				var columns = [];
//				columns.push(data[0]);
//				columns.push(data[1]);
//				for ( var i = 2; i < data.length; i++ )
//				{
//					if ( data[i] == 7 )	{	columns.push(true, true, true);		}
//					if ( data[i] == 6 )	{	columns.push(true, true, false);	}
//					if ( data[i] == 5 )	{	columns.push(true, false, true);	}
//					if ( data[i] == 4 )	{	columns.push(true, false, false);	}
//					if ( data[i] == 3 )	{	columns.push(false, true, true);	}
//					if ( data[i] == 2 )	{	columns.push(false, true, false);	}
//					if ( data[i] == 1 )	{	columns.push(false, false, true);	}
//					if ( data[i] == 0 )	{	columns.push(false, false, false);	}
//				}
//				
//
//				// SHOW COLUMNS: USERNAME, FULL NAME, PASSWORD
//				for ( var columnCounter = 0; columnCounter < columns.length; columnCounter++ )
//				{
//					var tableData = document.createElement("td");
//					
//					if ( columnCounter == 0 || columnCounter == 1 )
//					{
//						tableData.setAttribute('class', 'adminAccessTableData');	
//						var text = document.createTextNode(columns[columnCounter]);
//						tableData.appendChild(text);
//					}
//					else
//					{
//						if ( columns[columnCounter] == true )
//						{
//							tableData.setAttribute('class', 'adminAccessAllowed');	
//						}
//						else
//						{
//							tableData.setAttribute('class', 'adminAccessDenied');	
//						}
//						
//						// ADD ONCLICK TO TOGGLE CLASS
//						dojo.connect(tableData, "onclick", function(event){
//							adminAccess.togglePermission(event);
//						});
//					}
//
//					// APPEND TABLE DATA TO TABLE ROW
//					tableRow.appendChild(tableData);							
//				}
//
//				// APPEND TABLE DATA TO TABLE ROW
//				tableRow.appendChild(tableData);								
//				
//				// APPEND ROW TO TABLE
//				table.appendChild(tableRow);
//			}
//		}; // groupTable
//
//		this.groupStore.fetch( { query: {	project : "*" }, onComplete: groupTable});
//		
//	},	// buildTable
//
//
//	togglePermission : function (event)
//	{
//		
//		var target = event.target;
//		
//		var nodeClass = target.getAttribute('class');
//		if ( nodeClass == 'adminAccessAllowed' )
//		{
//			target.setAttribute('class', 'adminAccessDenied');
//		}
//		else
//		{
//			target.setAttribute('class', 'adminAccessAllowed');
//		}
//	},
//
//
//	//// MAKE SURE THAT PERMISSIONS ARE NUMERIC AND BETWEEN 1 AND 7
//	//checkPermissions : function (text)
//	//{
//	//	if ( ! text.match(/^\d+$/) ) return 1;
//	//	if ( text < 1 )	return 1;
//	//	if ( text > 7 ) return 1;
//	//	
//	//	return text;
//	//},
//
//	// CURRENTLY UNUSED CALLBACK FOR saveRow
//	updateAccess : function (json)
//	{
//	},
//	
//	
//	saveStore : function ()
//	{
//
//		// COLLECT DATA HERE
//		var dataArray = new Array;
//
//		var rights = [ "ownerrights", "grouprights", "worldrights" ];
//		var table = dojo.byId(this.tableId);
//		var rows = table.childNodes;
//		for ( var i = 2; i < rows.length; i++ )
//		{
//			var data = new Object;
//			var tableDatas = rows[i].childNodes;
//			data.project = tableDatas[0].innerHTML.toString();
//			var rightsCounter = 0;
//			for ( var j = 1; j < tableDatas.length; j++ )
//			{
//				var value = 0;
//				var values = [ 4, 2, 1];
//				for ( var k = 0; k < 3; k++ )
//				{
//					var node = tableDatas[(j + k)];
//					var nodeClass = node.getAttribute('class').toString();
//					if ( nodeClass == 'adminAccessAllowed' )
//					{
//						value += values[k];
//					}
//				}
//				
//				data[rights[rightsCounter]] = value;
//				rightsCounter++;
//				
//				j+=2;
//			}
//			
//			dataArray.push(data);
//		}
//
//		
//
//		// SAVE FOR LATER: DO IT USING THE STORE			
//		////////// print out for debugging
//		////////var adminAccess = this;
//		////////
//		////////var store = this.groupStore;
//		////////
//		////////for ( var i = 0; i < store._arrayOfTopLevelItems.length; i++ )
//		////////{
//		////////}
//		////////
//		////////// SEND TO SERVER AS RAW XHR POST
//		////////var dataArray = new Array;
//		////////for ( var i = 0; i < store._arrayOfTopLevelItems.length; i++ )
//		////////{
//		////////	
//		////////	var data = new Object;
//		////////	var item = store._arrayOfTopLevelItems[i];
//		////////	var attributes = store.getAttributes(item);
//		////////	if (attributes && attributes.length > 0)
//		////////	{
//		////////		for ( var j = 0; j < attributes.length; j++ )
//		////////		//for ( var j = 0; j < 1; j++ )
//		////////		{
//		////////			var values = store.getValues(item, attributes[j]);
//		////////
//		////////			if ( values )
//		////////			{
//		////////				
//		////////				// MULTI-VALUE ATTRIBUTE
//		////////				if (values.length > 1 )
//		////////				{
//		////////					data[attributes[i]] = [];
//		////////					for ( var k = 0; k < values.length; k++ )
//		////////					{
//		////////						var value = values[k];
//		////////						data[attributes[j]].push(value);
//		////////					}
//		////////				}
//		////////				// SINGLE VALUE ATTRIBUTE
//		////////				else
//		////////				{
//		////////					data[attributes[j]] = values[0];
//		////////				}
//		////////			}
//		////////		}
//		////////		
//		////////	}
//		////////	dataArray.push(data);
//		////////}
//
//		var url = Agua.cgiUrl + "agua?";
//
//		// CREATE JSON QUERY
//		var query = new Object;
//		query.username = Agua.cookie('username');
//		query.sessionId = Agua.cookie('sessionId');
//		query.mode = "saveAccess";
//		query.data = dojo.toJson(dataArray);
//		
//		// SEND TO SERVER
//		dojo.xhrPut(
//			{
//				url: url,
//				contentType: "text",
//				putData: dojo.toJson(query),
//				timeout: 3000,
//				load: function(response, ioArgs) {
//					return response;
//				},
//				error: function(response, ioArgs) {
//					return response;
//				}
//			}
//		);
//		
//	},
//
//	loadingFunction : function (widget, name, id){
//		
//	}
}

); // plugins.admin.Access
