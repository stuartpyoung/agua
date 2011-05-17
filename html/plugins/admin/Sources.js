dojo.provide("plugins.admin.Sources");

// ALLOW THE USER TO ADD, REMOVE AND MODIFY SOURCES

//dojo.require("dijit.dijit"); // optimize: load dijit layer
dojo.require("dijit.form.Button");
dojo.require("dijit.form.TextBox");
dojo.require("dijit.form.Textarea");
dojo.require("dojo.parser");
dojo.require("dojo.dnd.Source");
dojo.require("plugins.core.Common");

// HAS A
dojo.require("plugins.admin.SourceRow");

dojo.declare(
    "plugins.admin.Sources",
	[ dijit._Widget, dijit._Templated, plugins.core.Common ],
{
	//Path to the template of this widget. 
	templatePath: dojo.moduleUrl("plugins", "admin/templates/sources.html"),

	// Calls dijit._Templated.widgetsInTemplate
	widgetsInTemplate : true,

	//addingSource STATE
	addingSource : false,

	// OR USE @import IN HTML TEMPLATE
	cssFiles : [ "plugins/admin/css/sources.css" ],

	// PARENT WIDGET
	parentWidget : null,


	constructor : function(args)
	{
		// GET INFO FROM ARGS
		this.parentWidget = args.parentWidget;
		this.sources = args.parentWidget.sources;

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
		this.tabContainer.addChild(this.sourcesTab);
		this.tabContainer.selectChild(this.sourcesTab);

		// SET DRAG SOURCE - LIST OF SOURCES
		this.setTable();

		// SET NEW SOURCE FORM
		this.setNewSource();

		// SET TRASH DROP TARGET
		this.setTrash();	
	},



	// reloadSourceTabs
	//
	// RELOAD RELEVANT DISPLAYS IN GROUP-RELATED TABS
	//
	reloadSourceTabs : function ()
	{

		var tabPaneNames = ["plugins.admin.GroupSources"];
		for ( var i in tabPaneNames )
		{
			if ( this.parentWidget.paneWidgets[tabPaneNames[i]] != null )
			{
				this.parentWidget.paneWidgets[tabPaneNames[i]].reload();
			}
		}
	},


	clearValue : function (textarea)
	{

		if ( textarea.clicked == true ) return;

		textarea.clicked = true;
		textarea.value = '';
		textarea.focus();
	},

	setNewSource : function ()
	{

		// FOCUS ON WIDGET TO FIX ITS WIDTH
		this.newName.focus();

		// SET ADD SOURCE ONCLICK
		dojo.connect(this.addSourceButton, "onClick", dojo.hitch(this, "addSource"));

		// SET ONCLICK TO CANCEL DEFAULT TEXT
		dojo.connect(this.newName, "onclick", dojo.hitch(this, "clearValue", this.newName));
		dojo.connect(this.newName, "onfocus", dojo.hitch(this, "clearValue", this.newName));
		dojo.connect(this.newDescription, "onclick", dojo.hitch(this, "clearValue", this.newDescription));
		dojo.connect(this.newDescription, "onfocus", dojo.hitch(this, "clearValue", this.newDescription));
		dojo.connect(this.newLocation, "onclick", dojo.hitch(this, "clearValue", this.newLocation));
		dojo.connect(this.newLocation, "onfocus", dojo.hitch(this, "clearValue", this.newLocation));

	},



	setTable : function ()
	{

		// DELETE EXISTING TABLE CONTENT
		while ( this.sourcesTable.firstChild )
		{
			this.sourcesTable.removeChild(this.sourcesTable.firstChild);
		}

		var dataArray = new Array;
		var sourceArray = Agua.getSources();

		sourceArray = this.sortHasharray(sourceArray, 'name');

		// CHECK sourceArray IS NOT NULL OR EMPTY
		if ( sourceArray == null || sourceArray.length == 0 )
		{
			return;
		}

		// GENERATE dataArray TO INSERT INTO DND SOURCE TABLE
		for ( var j = 0; j < sourceArray.length; j++ )
		{
			var data = sourceArray[j];				
			data.toString = function () { return this.name; }
			dataArray.push( { data: data, type: ["draggableItem"] } );
		}

		// GENERATE DND SOURCE
		var dragSource = new dojo.dnd.Source(
			this.sourcesTable,
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
			if ( node.location == null )
			{
				node.location = '';
			}

			var sourceRow = new plugins.admin.SourceRow({
				name : node.name,
				description : node.description,
				location : node.location,
				parentWidget : this
			});

			node.innerHTML = '';
			node.appendChild(sourceRow.domNode);

		}

		var sourceObject = this;
		//dojo.connect(dragSource, "creator", sourceObject.specialAvatar );

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


	editSourceRow : function (sourceRow, node)
	{

		// RETURN IF ALREADY EDITING SOURCE ROW (I.E., MULTIPLE CLICKS)
		if ( this.editingSourceRow == true ) return;
		this.editingSourceRow = true;

		// REPLACE THE TD INNERHTML WITH A TEXTAREA
		var text = node.innerHTML;

		node.innerHTML = '';
		if ( text == null || ! text ) text = '';

		// RETURN IF THIS IS A DOUBLE-CLICK
		if ( text.match(/^<i/) ||
			text.match(/^<br/) ||
			text.match(/^<fieldset/) ||
			text.match(/^<textarea/) )
		{
			this.editingSourceRow = false;
			return;
		}

		// CREATE INPUT TEXT AREA
		var textarea = document.createElement('textarea');
		dojo.addClass(textarea, 'editSourceRow');
		node.appendChild(textarea);
		textarea.value = text;
		textarea.focus();

		// SET NEW PROJECT LISTENER
		var sourcesObject = this;
		dojo.connect(textarea, "onkeypress", function(evt){

			// summary: handles keyboard events
			var key = evt.charOrCode;

			if ( key == 13 )
			{
				var newText = textarea.value;

				var source = new Object;
				source.name = sourceRow.name.innerHTML;
				source.description = sourceRow.description.innerHTML;
				source.location = sourceRow.location.innerHTML;

				if ( source.description.match(/^<textarea/) )
					source.description = sourceRow.description.firstChild.value;

				if ( source.location.match(/^<textarea/) )
					source.location = sourceRow.location.firstChild.value;

				// REMOVE WHITESPACE
				source.name = source.name.replace(/^\s+/, '');
				source.name = source.name.replace(/\s+$/, '');
				source.description = source.description.replace(/^\s+/, '');
				source.description = source.description.replace(/\s+$/, '');
				source.location = source.location.replace(/^\s+/, '');
				source.location = source.location.replace(/\s+$/, '');


				if ( source.name != '' )
					//&& source.description != '' 
					//&& source.location != '' )
				{
					// REMOVE ORIGINAL source OBJECT FROM Agua.sources ARRAY
					// THEN ADD NEW SOURCE OBJECT TO Agua.sources ARRAY
					Agua.removeSource({ name: source.name});
					Agua.addSource(source);

					// REMOVE TEXTAREA
					node.removeChild(textarea);
					node.innerHTML = newText;

					// REDO SOURCE TABLE
					sourcesObject.setTable();

					// SAVE NEW SOURCE TO REMOTE DATABASE
					sourcesObject.saveSource(source);

					this.editingSourceRow = false;
				}
			}
		});

	},



	deleteSource : function (name)
	{

		// CLEAN UP WHITESPACE
		name = name.replace(/\s+$/,'');
		name = name.replace(/^\s+/,'');

		var sourceObject = { name: name };

		// REMOVING SOURCE FROM Agua.sources
		var success = Agua.removeSource(sourceObject)

		// RESET THE SOURCES TABLE
		this.setTable();

		var url = Agua.cgiUrl + "/agua?";

		// CREATE JSON QUERY
		var query = new Object;
		query.username = Agua.cookie('username');
		query.sessionId = Agua.cookie('sessionId');
		query.mode = "deleteSource";
		query.data = sourceObject;

		// SEND TO SERVER
		dojo.xhrPut(
			{
				url: url,
				contentType: "text",
				putData: dojo.toJson(query),
				timeout: 15000,
				load: function(response, ioArgs) {
					return response;
				},
				error: function(response, ioArgs) {
					return response;
				}
			}
		);


		// RELOAD SOURCE-RELATED TABS
		this.reloadSourceTabs();

	}, // Sources.deleteSource


	addSource : function (source)
	{

		var name = this.newName.value;
		var description = this.newDescription.value;
		var location = this.newLocation.value;

		// CHECK FOR VALID INPUTS
		if ( name == '' || name.match(/^\s*Name\s*$/) )
		{
			dojo.addClass(this.newName, 'invalid');
		}
		else{
			dojo.removeClass(this.newName, 'invalid');
		}
		if ( description == '' || description.match(/^\s*Description\s*$/) )
		{
			dojo.addClass(this.newDescription, 'invalid');
		}
		else{
			dojo.removeClass(this.newDescription, 'invalid');
		}
		if ( location == '' || location.match(/^\s*Location\s*$/) )
		{
			dojo.addClass(this.newLocation, 'invalid');
		}
		else{
			dojo.removeClass(this.newLocation, 'invalid');
		}

		if ( name == '' || name.match(/^\s*Name\s*$/) 
			|| description == '' || description.match(/^\s*Description\s*$/)
			|| location == '' || location.match(/^\s*Location\s*$/) )
		{
			return;
		}

		var sourceObject = { name: name, description: description, location: location };

		if ( Agua.isSource(sourceObject) )
		{
			return;
		}

		// ADD SOURCE TO SOURCES ARRAY
		Agua.addSource(sourceObject);

		// RESET SOURCE TABLE
		this.setTable();

		// SAVE SOURCE TO REMOTE DATABASE
		this.saveSource(sourceObject);

		// RELOAD SOURCE-RELATED TABS
		this.reloadSourceTabs();

	},	// Sources.addSource



	saveSource : function (source)
	{

		if ( this.savingSource == true )
		{
			return;
		}
		this.savingSource = true;

		// CLEAN UP WHITESPACE AND SUBSTITUTE NON-JSON SAFE CHARACTERS
		source.originalName = this.jsonSafe(source.originalName, 'toJson');
		source.name = this.jsonSafe(source.name, 'toJson');
		source.description = this.jsonSafe(source.description, 'toJson');
		source.location = this.jsonSafe(source.location, 'toJson');

		var url = Agua.cgiUrl + "/agua?";

		// CREATE JSON QUERY
		var query = new Object;
		query.username = Agua.cookie('username');
		query.sessionId = Agua.cookie('sessionId');
		query.mode = "saveSource";
		query.data = source;

		// SEND TO SERVER
		dojo.xhrPut(
			{
				url: url,
				contentType: "text",
				putData: dojo.toJson(query),
				timeout: 15000,
				load: function(response, ioArgs) {
					return response;
				},
				error: function(response, ioArgs) {
					return response;
				}
			}
		);

		this.savingSource = false;

	}, // Sources.saveSource




	// setTrash
	//
	//	DELETE NODE IF DROPPED INTO TRASH. ACTUAL REMOVAL FROM THE
	//	DATA IS ACCOMPLISHED IN THE onDndDrop LISTENER OF THE SOURCE
	//
	setTrash : function ()
	{

		var trash = new dojo.dnd.Source(
			this.trashContainer,
			{
				accept : [ "draggableItem" ]
			}
		);

		// REMOVE DUPLICATE NODES
		var thisObject = this;
		dojo.connect(trash, "onDndDrop", function(source, nodes, copy, target){
			// NODE DROPPED ON SELF --> DELETE THE NODE
			if ( target == this )
			{

				for ( var i = 0; i < nodes.length; i++ )
				{
					var node = nodes[i];

					// HACK TO AVOID THIS ERROR: node.parentNode is null
					try {
						node.parentNode.removeChild(node);

						thisObject.deleteSource(node.name);
					}
					catch (e) {
					}
				}

				// DELETE EXISTING TABLE CONTENT
				//console.dir(thisObject.trashContainer);
				while ( thisObject.trashContainer.childNodes.length > 2 )
				{

					thisObject.trashContainer.removeChild(thisObject.trashContainer.childNodes[2]);
				}


			}
		});
	}

}); // plugins.admin.Sources

