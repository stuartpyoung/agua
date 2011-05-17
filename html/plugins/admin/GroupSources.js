dojo.provide("plugins.admin.GroupSources");

dojo.require("dijit.dijit"); // optimize: load dijit layer

dojo.require("dijit.form.Button");
dojo.require("dijit.form.ComboBox");
dojo.require("plugins.core.ComboBox");

dojo.require("dojo.data.ItemFileWriteStore");

// TOOLTIP
dojo.require("dijit.Tooltip");

dojo.require("dojo.parser");
dojo.require("dojo.dnd.Source");

// HAS A
dojo.require("plugins.admin.GroupSourceRow");

dojo.declare(
    "plugins.admin.GroupSources",

	[ dijit._Widget, dijit._Templated ],
{
	//Path to the template of this widget. 
	templatePath: dojo.moduleUrl("plugins", "admin/templates/groupsources.html"),

	// Calls dijit._Templated.widgetsInTemplate
	widgetsInTemplate : true,

	//addingSource STATE
	addingSource : false,

	// OR USE @import IN HTML TEMPLATE
	cssFiles : [ "plugins/admin/css/groupsources.css" ],

	// PARENT WIDGET
	parentWidget : null,

 	// CONSTRUCTOR	
	constructor : function(args) {

		// REGISTER module path FOR PLUGINS


        // LOAD SORIA AND FILEPICKER CSS
        this.loadCSS();		
	},

	postCreate : function()
	{

		this.startup();
	},


	startup : function ()
	{

		// COMPLETE CONSTRUCTION OF OBJECT
		this.inherited(arguments);	 

		// ADD ADMIN TAB TO TAB CONTAINER		
		this.parentWidget.rightTabContainer.addChild(this.sourcesTab);
		this.parentWidget.rightTabContainer.selectChild(this.sourcesTab);

		// SET GROUP COMBO
		this.setGroupCombo();

		// SET DRAG SOURCE - LIST OF SOURCES
		this.setDragSource();

		// SET DROP TARGET - SOURCES ALREADY IN THE GROUP
		this.setDropTarget();

		// SET TRASH DROP TARGET
		this.setTrash();


		//// example subscribe to events
		//dojo.subscribe("/dnd/start", function(source){
		//	//console.debug("/dnd/start Starting the drop", source);
		//});
		//dojo.subscribe("/dnd/drop/before", function(source, nodes, copy, target){
		//	//console.debug(copy ? "/dnd/drop/before Copying from" : "Moving from", source, "to", target, "before", target.before);
		//});
		//dojo.subscribe("/dnd/drop", function(source, nodes, copy, target){
		//	//console.debug(copy ? "/dnd/drop Copying from" : "Moving from", source, "to", target, "before", target.before);
		//});

	},


	// reload
	//
	// RELOAD THE GROUP COMBO AND DRAG SOURCE
	// (CALLED AFTER CHANGES TO SOURCES OR GROUPS DATA IN OTHER TABS)
	//
	reload : function ()
	{
		// SET GROUP COMBO
		this.setGroupCombo();

		// SET DRAG SOURCE - LIST OF SOURCES
		this.setDragSource();

		// SET DROP TARGET - SOURCES ALREADY IN THE GROUP
		this.setDropTarget();
	},


	// SET GROUP COMBO BOX
	setGroupCombo : function ()
	{

		// GET GROUP NAMES		
		var groupNames = Agua.getGroupNames();

		// SET STORE
		var data = {identifier: "name", items: []};
		for ( var i = 0; i < groupNames.length; i++ )
		{
			data.items[i] = { name: groupNames[i]	};
		}
		var store = new dojo.data.ItemFileWriteStore(	{	data: data	}	);

		this.groupCombo.popupClass = "groupsources groupCombo dijitReset dijitMenu";
		this.groupCombo.wrapperClass = "groupsources dijitPopup";
		this.groupCombo.itemHeight = 30;


		// SET COMBO
		this.groupCombo.store = store;
		this.groupCombo.startup();

		// SET COMBO VALUE
		var firstValue = groupNames[0];
		this.groupCombo.setValue(firstValue);

		// CONNECT ONCLICK WITH dojo.connect TO BUILD TABLE
		var adminSources = this;
		dojo.connect(this.groupCombo, "onChange", function(event) {
			adminSources.setDropTarget();
		});
	},



	setDragSource : function ()
	{
		// REMOVE ALL EXISTING CONTENT
		while ( this.dragSourceContainer.firstChild )
		{
			if ( dijit.byNode(this.dragSourceContainer)
				&& dijit.byNode(this.dragSourceContainer).destroy )
			{
				dijit.byNode(this.dragSourceContainer).destroy();
			}
			this.dragSourceContainer.removeChild(this.dragSourceContainer.firstChild);
		}

		var dataArray = new Array;
		var sourceArray = Agua.getSources();


		if ( sourceArray == null || ! sourceArray )
		{
			return;
		}

		for ( var j = 0; j < sourceArray.length; j++ )
		{
			var data = sourceArray[j];				
			data.toString = function () { return this.name; }
			dataArray.push( { data: data, type: ["draggableItem"] } );
		}

		// GENERATE DND SOURCE
		var dragSource = new dojo.dnd.Source(
			this.dragSourceContainer,
			{
				copyOnly: true,
				selfAccept: false,
				accept : [ "none" ]
			}
		);
		dragSource.insertNodes(false, dataArray);

		// SET TABLE ROW STYLE IN dojDndItems
		var allNodes = dragSource.getAllNodes();
		for ( var k = 0; k < allNodes.length; k++ )
		{
			// ADD CLASS FROM type TO NODE
			var node = allNodes[k];


			// SET NODE name AND description
			node.name = dataArray[k].data.name;
			node.description = dataArray[k].data.description;
			node.location = dataArray[k].data.location;

			var source = {
				name : node.name,
				location : node.location,
				description : node.description
			};

			source.parentWidget = this;			

			var groupSourceRow = new plugins.admin.GroupSourceRow(source);
			node.innerHTML = '';
			node.appendChild(groupSourceRow.domNode);
		}

		var groupObject = this;
		dragSource.creator = function (item, hint)
		{

			var node = dojo.doc.createElement("div");
			node.name = item.name;
			node.description = item.description;
			node.location = item.location;
			node.id = dojo.dnd.getUniqueId();
			node.className = "dojoDndItem";


			// SET FANCY FORMAT IN NODE INNERHTML
			node.innerHTML = "<table> <tr><td><strong style='color: darkred'>" + item.name + "</strong></td></tr><tr><td> " + item.description + "</td></tr></table>";

			return {node: node, data: item, type: ["text"]};
		};
	},


	setDropTarget : function ()
	{

		// DELETE EXISTING CONTENTS OF DROP TARGET
		while ( this.dropTargetContainer.firstChild )
		{
			this.dropTargetContainer.removeChild(this.dropTargetContainer.firstChild);
		}

		//// GET THE SOURCES IN THIS GROUP
		var groupname = this.groupCombo.getValue();
		var sourceArray = Agua.getSourcesByGroup(groupname);
		if ( sourceArray == null )
		{
			return;
		}		

		var dataArray = new Array;
		for ( var j = 0; j < sourceArray.length; j++ )
		{
			var data = sourceArray[j];				
			data.toString = function () { return this.name; }
			var description = sourceArray[j].description;

			dataArray.push( { data: data, type: ["draggableItem"], description: description  } );
		}

		// GENERATE DROP TARGET
		var dropTarget = new dojo.dnd.Source(
			this.dropTargetContainer,
			{
				//copyOnly: true,
				//selfAccept: false
				//,
				accept : [ "draggableItem" ]
			}
		);
		dropTarget.insertNodes(false, dataArray );

		// SET TABLE ROW STYLE IN dojDndItems
		var allNodes = dropTarget.getAllNodes();
		for ( var k = 0; k < allNodes.length; k++ )
		{
			// ADD CLASS FROM type TO NODE
			var node = allNodes[k];
			var nodeClass = dataArray[k].type;
			dojo.addClass(node, nodeClass);



			// SET APPLICATION FOR NODE
			node.description = dataArray[k].description;
			node.location = dataArray[k].location;

			var name = dataArray[k].data;
			node.name = name;


			// SET NODE name AND description
			node.name = dataArray[k].data.name;
			node.description = dataArray[k].data.description;
			node.location = dataArray[k].data.location;

			var source = {
				name : node.name,
				location : node.location,
				description : node.description
			};
			source.parentWidget = this;			

			var groupSourceRow = new plugins.admin.GroupSourceRow(source);
			node.innerHTML = '';
			node.appendChild(groupSourceRow.domNode);
		}

		var adminSource = this;
		dropTarget.creator = function (item, hint)
		{

			var name = item.name;
			var description = item.description;

			// ADD VALUES TO NODE SO THAT THEY GET PASSED TO this.addToGroup
			var node = dojo.doc.createElement("div");
			node.name = item.name;
			node.description = item.description;
			node.location = item.location;
			node.id = dojo.dnd.getUniqueId();
			node.className = "dojoDndItem";

			// SET FANCY FORMAT IN NODE INNERHTML
			//node.innerHTML = "<table> <tr><td><strong style='color: darkred'>" + name + "</strong></td></tr><tr><td> " + description + "</td></tr></table>";

			var source = {
				name : node.name,
				description : node.description,
				location : node.location
			};

			source.parentWidget = this;			

			var groupSourceRow = new plugins.admin.GroupSourceRow(source);
			node.innerHTML = '';
			node.appendChild(groupSourceRow.domNode);


			if ( hint != 'avatar' )
			{
				adminSource.addToGroup(node.name, node.description, node.location);
			}

			return {node: node, data: item, type: ["draggableItem"]};
		};


		dojo.connect(dropTarget, "onDndDrop", function(source, nodes, copy, target){

			if ( source == this )
			{
				for ( var i = 0; i < nodes.length; i++ )
				{
					var node = nodes[i];

					adminSource.removeFromGroup(node.name, node.description, node.location);
				}
			}

			// REMOVE DUPLICATE NODES
			var currentNodes = dropTarget.getAllNodes();
			if ( currentNodes == null || ! currentNodes )
			{
				return;
			}

			var names = new Object;
			for ( var i = 0; i < currentNodes.length; i++ )
			{
				var node = currentNodes[i];
				if ( ! names[node.name] )
				{
					names[node.name] = 1;
				}
				else
				{
					// HACK TO AVOID THIS ERROR: node.parentNode is null
					try {
						node.parentNode.removeChild(node)

						// NB: this has no effect on nodes in dnd.Source
						//currentNodes.splice(i, 1); 
						//i--;
					}
					catch (e) {
					}
				}
			}
		});
	},

	addToGroup : function (name, description, location )
	{

		var groupObject = new Object;
		groupObject.username = Agua.cookie('username');
		groupObject.name = name;
		groupObject.description = description;
		groupObject.location = location;

		// ADD SOURCE OBJECT TO THE SOURCES IN THIS GROUP
		var groupName = this.groupCombo.getValue();

		if ( Agua.addSourceToGroup(groupName, groupObject) == false )
		{
			return;
		}

		// ADD THE SOURCE INTO THE groupusers TABLE ON THE SERVER
		var data = new Object;
		data.name = name;
		data.description = description;
		data.location = location;
		data.groupname = groupName;
		data.type = "source";

		var url = Agua.cgiUrl + "agua?";
		var query = new Object;
		query.username = Agua.cookie('username');
		query.sessionId = Agua.cookie('sessionId');
		query.data = data;
		query.mode = "addToGroup";

		var queryString = dojo.toJson(query).replace(/undefined/g, '""');

		// SEND TO SERVER
		dojo.xhrPut(
			{
				url: url,
				contentType: "text",
				sync : false,
				handleAs: "json",
				putData: queryString,
				timeout: 15000,
				load: function(data)
				{
				},
				error: function(response, ioArgs) {
					return response;
				}
			}
		);
	},



	removeFromGroup : function (name, description, location )
	{

		// REMOVE SOURCE FROM THE SOURCES IN THIS GROUP
		var groupName = this.groupCombo.getValue();

		var groupObject = new Object;
		groupObject.name = name;
		groupObject.description = description;
		groupObject.location = location;

		if ( Agua.removeSourceFromGroup(groupName, groupObject) == false )
		{
			return;
		}

		// REMOVE THE SOURCE FROM THE groupusers TABLE ON THE SERVER
		var data = new Object;
		data.name = String(name);
		data.description = description;
		data.location = location;
		data.groupname = groupName;
		data.type = "source";

		var url = Agua.cgiUrl + "agua?";
		var query = new Object;
		query.username = Agua.cookie('username');
		query.sessionId = Agua.cookie('sessionId');
		query.data = data;
		query.mode = "removeFromGroup";

		var queryString = dojo.toJson(query).replace(/undefined/g, '""');

		// SEND TO SERVER
		dojo.xhrPut(
			{
				url: url,
				contentType: "text",
				sync : false,
				handleAs: "json",
				putData: queryString,
				timeout: 15000,
				load: function(data)
				{
				},
				error: function(response, ioArgs) {
					return response;
				}
			}
		);
	},


    loadCSS : function ()
    {

		// LOAD CSS
		for ( var i in this.cssFiles )
		{

			var cssFile = this.cssFiles[i];
			var cssNode = document.createElement('link');
			cssNode.type = 'text/css';
			cssNode.rel = 'stylesheet';
			cssNode.href = cssFile;
		cssNode.id = "themeStyles";
			document.getElementsByTagName("head")[0].appendChild(cssNode);
		}

    },


	setTrash : function ()
	{

		var trash = new dojo.dnd.Source(
			this.trashContainer,
			{
				accept : [ "draggableItem" ]
			}
		);

		// REMOVE DUPLICATE NODES
		dojo.connect(trash, "onDndDrop", function(source, nodes, copy, target){


			// NODE DROPPED ON SELF --> DELETE THE NODE
			if ( target == this )
			{
				var currentNodes = trash.getAllNodes();

				for ( var i = 0; i < currentNodes.length; i++ )
				{
					var node = currentNodes[i];


					// HACK TO AVOID THIS ERROR: node.parentNode is null
					try {
						node.parentNode.removeChild(node)
					}
					catch (e) {
					}
				}
			}
		});
	}


}

); // plugins.admin.GroupSources

