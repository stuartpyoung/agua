dojo.provide("plugins.workflow.SharedApps");

// NB: ADMIN USER MAY ADD, REMOVE AND MODIFY APPS IN 'Applications'
// NB: ORDINARY USER MAY DO THE SAME FOR APPS IN 'Custom'

//dojo.require("dijit.dijit"); // optimize: load dijit layer
dojo.require("dijit.form.Button");
dojo.require("dijit.form.TextBox");
dojo.require("dijit.form.Textarea");
dojo.require("dojo.parser");
dojo.require("dojo.dnd.Source");
dojo.require("plugins.core.Common");

// HAS A
dojo.require("plugins.workflow.AppType");
dojo.require("plugins.workflow.AppRow");
dojo.require("plugins.workflow.AppsMenu");

dojo.declare(
    "plugins.workflow.SharedApps",
	[ dijit._Widget, dijit._Templated, plugins.core.Common ],
{

//Path to the template of this widget. 
templatePath: dojo.moduleUrl("plugins", "workflow/templates/sharedapps.html"),

// Calls dijit._Templated.widgetsInTemplate
widgetsInTemplate : true,

//addingApp STATE
addingApp : false,

// OR USE @import IN HTML TEMPLATE
//cssFiles : [ dojo.moduleUrl("plugins") + "/workflow/css/apps.css" ],
cssFiles : [ dojo.moduleUrl("plugins") + "/workflow/css/apps.css"],

// PARENT WIDGET
parentWidget : null,

// TAB CONTAINER
tabContainer : null,

// CONTEXT MENU
contextMenu : null,

// CORE WORKFLOW OBJECTS
core : null,

/////}

constructor : function(args) {
	this.core = args.core;

	// GET INFO FROM ARGS
	this.parentWidget = args.parentWidget;
	this.tabContainer = args.tabContainer;

	// LOAD SORIA AND FILEPICKER CSS
	this.loadCSS();		
},

postCreate : function() {

	this.startup();
},


startup : function () {

	// COMPLETE CONSTRUCTION OF OBJECT
	this.inherited(arguments);	 

	// ADD TO TAB CONTAINER		

	this.tabContainer.addChild(this.mainTab);
	this.tabContainer.selectChild(this.mainTab);

	// CREATE SOURCE MENU
	this.setContextMenu();

	// SET DRAG APP - LIST OF APPS
	this.setDragSource();

	// SUBSCRIBE TO UPDATES
	Agua.updater.subscribe(this, "updateApps");
},

// GENERATE CONTEXT MENU
setContextMenu : function () {

	this.contextMenu = new plugins.workflow.AppsMenu(
		{
			parentWidget: this
		}
	);
},

// REFRESH APPS
refresh : function () {

	this.updateApps();
},

updateApps : function (args) {

	this.setDragSource();
},

setDragSource : function () {

	// DELETE EXISTING TABLE CONTENT
	while ( this.appsNode.firstChild )
	{
		this.appsNode.removeChild(this.appsNode.firstChild);
	}
	if ( Agua.hasApps("common") == false )	return;

	var types = Agua.getAppTypes("common");

	for ( var i = 0; i < types.length; i++ )
	{
		// GET APPLICATIONS
		var sourceArray = Agua.getAppsByType(types[i], "common");

		// CHECK sourceArray IS NOT NULL OR EMPTY
		if ( sourceArray == null || sourceArray.length == 0 )
		{
			return;
		}

		// CREATE TITLE PANE
		var titlePane = new plugins.workflow.AppType(
		{
			title: types[i]
		});
		this.appsNode.appendChild(titlePane.domNode);


		// GENERATE dataArray TO INSERT INTO DND APP TABLE
		var dataArray = new Array;
		for ( var j = 0; j < sourceArray.length; j++ )
		{
			var data = sourceArray[j];				
			dataArray.push( { data: data, type: ["draggableItem"] } );
		}

		// GENERATE DND APP
		var dragSource = new dojo.dnd.Source(
			titlePane.dragSource,
			{
				copyOnly: true,
				selfAccept: false,
				accept : [ "none" ]
			}
		);
		dragSource.insertNodes(false, dataArray);

		// SET TABLE ROW STYLE IN dojDndItems
		var allNodes = dragSource.getAllNodes();
		for ( var k = 0; k < sourceArray.length; k++ )
		{
			// ADD CLASS FROM type TO NODE
			var node = allNodes[k];


			// SET NODE name AND description
			node.name = dataArray[k].data.name;
			node.owner = dataArray[k].data.owner;
			node.type = dataArray[k].data.type;
			node.executor = dataArray[k].data.executor;
			node.localonly = dataArray[k].data.localonly;
			node.location = dataArray[k].data.location;
			node.description = dataArray[k].data.description;
			node.notes = dataArray[k].data.notes;
			if ( node.description == null || node.description == '' )
			{
				node.description = "Description";
			}
			if ( node.notes == null || node.notes == '' )
			{
				node.notes = "Notes";
			}

			// SET SOURCE HASH USED TO INSTANTIATE plugins.workflow.AppRow
			var source = {
				name : node.name,
				owner : node.owner,
				type : node.type,
				executor : node.executor,
				localonly : node.localonly,
				location : node.location,
				description : node.description,
				notes : node.notes
			};
			// ADD TO node.application FOR INTERFACE WITH Workflow.setDropTarget
			node.application = source;

			//source.parentWidget = this;

			// ADD SOURCE MENU
			this.contextMenu.bind(node);

			// INSTANTIATE plugins.workflow.AppRow AND APPEND TO NODE
			var sourceRow = new plugins.workflow.AppRow(source);


			node.innerHTML = '';
			node.appendChild(sourceRow.domNode);
		}

		var sourceObject = this;
		dragSource.creator = function (item, hint)
		{

			var node = dojo.doc.createElement("div");
			node.application = item;
			node.application.dummy = "dummy";
			node.name = item.name;
			node.owner = item.owner;
			node.description = item.description;
			node.notes = item.notes;
			node.id = dojo.dnd.getUniqueId();
			node.className = "dojoDndItem";


			// SET FANCY FORMAT IN NODE INNERHTML
			node.innerHTML = "<table> <tr><td><strong style='color: darkred'>" + item.name + "</strong></td></tr><tr><td> " + item.description + "</td></tr></table>";


			// SET SOURCE HASH USED TO INSTANTIATE plugins.workflow.AppRow
			var source = {
				name : item.name,
				owner : item.owner,
				type : item.type,
				executor : item.executor,
				localonly : item.localonly,
				location : item.location,
				description : item.description,
				notes : item.notes
			};
			// ADD TO node.application FOR INTERFACE WITH Workflow.setDropTarget
			node.application = source;
			source.parentWidget = sourceObject;


			// INSTANTIATE plugins.workflow.AppRow AND APPEND TO NODE
			var sourceRow = new plugins.workflow.AppRow(source);
			node.innerHTML = '';
			node.appendChild(sourceRow.domNode);

			node.parentWidget = sourceRow;

			return {node: node, data: item, type: ["text"]};
		};
	}
}



}); // plugins.workflow.SharedApps
