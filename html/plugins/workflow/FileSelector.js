dojo.provide("plugins.workflow.FileSelector");

dojo.require("dojox.widget.RollingList");

dojo.require("plugins.workflow._GroupSelectorPane");


dojo.declare("plugins.workflow._FileInfoPane", 
	 dojox.widget._RollingListPane, {
	// summary: a pane to display the information for the currently-selected
	//	file

	// templateString: string
	//	delete our template string
	//templateString: "",

	// templatePath: string
	//	Our template path
	templateString:"<div class=\"dojoxFileInfoPane\">\r\n\t<table>\r\n\t\t<tbody>\r\n\t\t\t<tr>\r\n\t\t\t\t<td class=\"dojoxFileInfoLabel dojoxFileInfoNameLabel\">Filename </td>\r\n\t\t\t\t<td class=\"dojoxFileInfoName\" dojoAttachPoint=\"nameNode\"></td>\r\n\t\t\t</tr>\r\n\t\t\t<tr>\r\n\t\t\t\t<td class=\"dojoxFileInfoLabel dojoxFileInfoSizeLabel\">Size (bytes) </td>\r\n\t\t\t\t<td class=\"dojoxFileInfoSize\" dojoAttachPoint=\"sizeNode\"></td>\r\n\t\t\t</tr>\r\n\t\t\t<tr>\r\n\t\t\t\t<td class=\"dojoxFileInfoLabel dojoxFileInfoPathLabel\">Sample </td>\r\n\t\t\t\t<td class=\"dojoxFileInfoPath\" dojoAttachPoint=\"sampleNode\"></td>\r\n\t\t\t</tr>\r\n\t\t</tbody>\r\n\t</table>\r\n</div>\r\n",

	postMixInProperties: function(){
		//this._messages = dojo.i18n.getLocalization("dojox.widget", "FileSelector", this.lang);
		//this.inherited(arguments);
	},

	onItems: function()
	{
		// summary:
		//	called after a fetch or load - at this point, this.items should be
		//  set and loaded.
		var store = this.store, item = this.items[0];




		if(!item){
			this._onError("Load", new Error("No item defined"));
		}else{
			this.nameNode.innerHTML = store.getLabel(item);
			this.sampleNode.innerHTML = store.getValue(item, "sample");
			this.sizeNode.innerHTML = store.getValue(item, "size");
			this.parentWidget.scrollIntoView(this);
			this.inherited(arguments);
		}
	}
});


dojo.declare("plugins.workflow.FileSelector", dojox.widget.RollingList, {
	// summary: a specialized version of RollingList that handles file information
	//  in a store

	className: "dojoxFileSelector",

	// pathSeparator: string
	//  Our file separator - it will be guessed if not set
	pathSeparator: "",

	// topDir: string
	//	The top directory string - it will be guessed if not set
	topDir: "",

	// parentAttr: string
	//	the attribute to read for finding our parent directory
	parentAttr: "parentDir",

	// pathAttr: string
	//  the attribute to read for getting the full path of our file
	pathAttr: "path",

	// PARENT WIDGET, I.E., FileManager
	parentWidget : null,

	// WORKFLOW WIDGET
	workflowObject : null,

	_itemsMatch: function(/*item*/ item1, /*item*/ item2){
		// Summary: returns whether or not the two items match - checks ID if
		//  they aren't the exact same object - ignoring trailing slashes


		return true;

	},


	// callback FROM OBJECT THAT GENERATED THE FileManager WHICH CREATED THIS FileSelector
	callback : null,


	// SET callback IN CONSTRUCTOR
	constructor : function(args) {


		this.callback = args.callback;
		this.parentWidget = args.parentWidget;
		this.workflowObject = args.workflowObject;


		//this.inherited(arguments);
	},



	startup: function(){


		if (this._started )
		{
			return;
		}

		this.inherited(arguments);

		// Figure out our file separator if we don't have it yet
		var conn, child = this.getChildren()[0];



		var setSeparator = dojo.hitch(this, function(){

			if(conn){
				this.disconnect(conn);
			}
			delete conn;

			var item = child.items[0];



			if(item){
				var store = this.store;


				var parent = store.getValue(item, this.parentAttr);
				var path = store.getValue(item, this.pathAttr);


				this.pathSeparator = this.pathSeparator || store.pathSeparator;
				if(!this.pathSeparator){
					this.pathSeparator = path.substring(parent.length, parent.length + 1);
				}
				if(!this.topDir){
					this.topDir = parent;
					if(this.topDir.lastIndexOf(this.pathSeparator) != (this.topDir.length - 1)){
						this.topDir += this.pathSeparator;
					}
				}


			}
		});
		if(!this.pathSeparator || !this.topDir){
			if(!child.items){
				conn = this.connect(child, "onItems", setSeparator);
			}else{
				setSeparator();
			}
		}
	},


	_removeAfter: function(/*Widget or int*/ idx){

		// summary: removes all widgets after the given widget (or index)
		if ( idx < 0 ) {	return;	}

		if(typeof idx != "number"){
			idx = this.getIndexOfChild(idx);
		}

		var numberChildren = this.getChildren().length;

			dojo.forEach(this.getChildren(), function(c, i){

				if(i > idx )
				{
					this.removeChild(c);
					c.destroy();
				}
			}, this);
		//}
		var children = this.getChildren(), child = children[children.length - 1];
		var selItem = null;
		while(child && !selItem){
			var val = child._getSelected ? child._getSelected() : null;
			if(val){
				selItem = val.item;
			}
			child = child.parentPane;
		}
		if(!this._setInProgress){
			this._setValue(selItem);
		}
	},


	addChild: function(/*Widget*/ widget, /*int?*/ insertIndex){


		// summary: adds a child to this rolling list - if passed an insertIndex,
		//  then all children from that index on will be removed and destroyed
		//  before adding the child.


		//if ( insertIndex > 0 )
		//{
			this._removeAfter(insertIndex - 1);
		//}
		//else
		//{
		//	this._removeAfter(0);
		//}

		this.inherited(arguments);



		if(!widget._started){
			widget.startup();
			}

		this.layout();

		if(!this._savedFocus){
			widget.focus();
		}
	},




	getChildItems: function(item){


		var ret = this.inherited(arguments);


		// CHECK IF THE ITEM IS A DIRECTORY AND EMPTY
		if(!ret && this.store.getValue(item, "directory")){
			// It's an empty directory - so pass through an empty array
			ret = [];
		}

		return ret;
	},


	_getMenuItemForItem: function(/*item*/ item, /* dijit._Contained */ parentPane)
	{
	// summary: returns a widget for the given store item.  The returned
	//  item will be added to this widget's container widget.  null will
	//  be passed in for an "empty" item.


		if ( item )
		{
			if ( item.children )
			{
				for ( var i = 0; i < item.children.length; i++ )
				{
				}
			}
		}

			var widgetItem;
			widgetItem = this.getMenuItemForItem(item, parentPane, null);

			// DO THIS LATER
			//this._updateClass(widgetItem.domNode, "Item", {"Single": true});		

			// ADD this.store AND item TO WIDGET ITEM
			widgetItem.store = this.store;

			widgetItem.item = item;


			return widgetItem;
		//}
	},

	getMenuItemForItem: function(/*item*/ item, /* dijit._Contained */ parentPane, /* item[]? */ children){

		if ( ! item )
		{
			var menuItem = {
				data : 'EMPTY',
				type : [ 'file' ]
			};
			return menuItem;	
		}

		var parentWidgetNode = parentPane.parentWidget.domNode;
		var childNodes = parentWidgetNode.childNodes;

		var dragType = "workflow";
		if ( childNodes.length > 1 )
		{
			dragType = "file";
		}

		//var iconClass = "dojoxDirectoryItemIcon";
		//if(!this.store.getValue(item, "directory")){
		//	iconClass = "dojoxFileItemIcon";
		//	var l = this.store.getLabel(item), idx = l.lastIndexOf(".");
		//	if(idx >= 0){
		//		iconClass += " dojoxFileItemIcon_" + l.substring(idx + 1);
		//	}
		//}

		// DO LOGIC TO 
		// RETURN A DND SOURCE ITEM
		var menuItem = {
			data : item.name,
			type : [ dragType ]
		};


		return menuItem;
	},

	// OVERRIDE OF RollingList METHOD
	_getPaneForItem: function(/* item? */ item, /* dijit._Contained? */ parentPane, /* item[]? */ children){		// summary: gets the pane for the given item, and mixes in our needed parts
		// Returns the pane for the given item (null if the root pane) - after mixing in
		// its stuff.

		if ( item )
		{
		}
		else
		{
		}


		var pane = this.getPaneForItem(item, parentPane, children);

		// REQUIRED: SET 'store' OF RETURNED PANE FOR ITEM
		pane.store = this.store;


		this.store.path = this.path;
		if ( item && item.parentPath)
		{


			// NEEDED TO ADD CURRENT DIRECTORY TO PATH OF pane.store.path
			// I.E., pane.store.path = Project1/Workflow1/inputdir INSTEAD OF JUST Project1/inputdir
			// IN FETCH
			pane.store.path = item.parentPath;

			// NEEDED TO UPDATE FULL PATH IN FileStore's this.query
			// OTHERWISE WILL SUCCESSFULLY LOAD THE SAME BASE DIRECTORY (E.G., Project1) AGAIN AND AGAIN
			this.store.path += "/" + item.name

			pane.path = item.parentPath + "/" + item.name;
			//this.path + "/" + item.name;
		}
		else
		{
			pane.path = this.path;
		}



		pane.store.parentPath = this.path;

		pane.parentWidget = this;
		pane.parentPane = parentPane||null;

		if( !item )
		{
			pane.query = this.query;
			pane.queryOptions = this.queryOptions;
		}
		else if ( children )
		{
			pane.items = children;
		}
		else
		{
			pane.items = [item];
		}


		return pane;
	},



	getPaneForItem: function(/*item*/ item, /* dijit._Contained */ parentPane, /* item[]? */ children)
	{


		//if ( item )
		//{
		//}


		var groupDragPane = null;

		if(!item || (this.store.isItem(item) && children) || this.store.getValue(item, "directory") )
		//if(!item || (this.store.isItem(item) && this.store.getValue(item, "directory")))
		//if( ! item || item.children )
		{

			var path = '';;
			var parentPath = '';
			if ( item )
			{
				path = item.path;
				parentPath = item.parentPath;
			}



			groupDragPane = new plugins.workflow._GroupSelectorPane(
				{
					parentPath: parentPath,
					path: path,
					callback: this.callback,
					parentWidget : this,
					grandParentWidget : this.parentWidget,
					workflowObject: this.workflowObject
				}
			);



		}
		//else if(this.store.isItem(item) && !this.store.getValue(item, "directory"))
		else if(this.store.isItem(item) && !this.store.getValue(item, "directory"))
		{

			// GENERATE A FILE INFO PANE
			groupDragPane = new plugins.workflow._FileInfoPane({});
		}


		return groupDragPane;
	},

	_setPathValueAttr: function(/*string*/ path){
		// Summary: sets the value of this widget based off the given path
		if(!path){
			this.attr("value", null);
			return;
		}
		if(path.lastIndexOf(this.pathSeparator) == (path.length - 1)){
			path = path.substring(0, path.length - 1);
		}
		this.store.fetchItemByIdentity({identity: path,
										onItem: dojo.hitch(this, "attr", "value"),
										scope: this});
	},

	_getPathValueAttr: function(/*item?*/val){
		// summary: returns the path value of the given value (or current value
		//  if not passed a value)
		if(!val){
			val = this.value;
		}
		if(val && this.store.isItem(val)){
			return this.store.getValue(val, this.pathAttr);
		}else{
			return "";
		}
	}
});
