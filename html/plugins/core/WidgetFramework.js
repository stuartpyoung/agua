dojo.provide( "plugins.core.WidgetFramework");

// OBJECT:  WidgetFramework
// PURPOSE: CREATE A WIDGET FRAMEWORK

// LAYOUT CLASSES
dojo.require("dijit.layout.BorderContainer");
dojo.require("dijit.layout.LayoutContainer");
dojo.require("dijit.layout.AccordionContainer");
dojo.require("dijit.layout.ContentPane");
dojo.require("dijit.layout.SplitContainer");
dojo.require("dijit.layout.TabContainer");
dojo.require("dijit.form.Button");
dojo.require("dijit.Toolbar");
dojo.require("dojox.layout.ExpandoPane");

dojo.declare( "plugins.core.WidgetFramework", null,
{
	// widget counter   
	cnt: 1,

	// create the widgets specified in layout and add them to widget "rootWidget"
	create : function (layout, args)
	{
		var rootWidgetId = args.root;
		var deleteRoot = args.deleteRoot;

		if ( deleteRoot != false && ! deleteRoot )
		//if ( deleteRoot == null || ! deleteRoot )
		{
			deleteRoot = true;
		}

		// erase old widget hierarchy (if it exists and deleteRoot == true)
		var rootWidget = dijit.byId(rootWidgetId);
		if(rootWidget && deleteRoot)
		{
			rootWidget.destroyRecursive();
		}

		// create new widget
		rootWidget = this.createWidgetHierarchy(layout);

		// and display it
		var wrapper = dojo.byId(args.target);

		if ( wrapper == null || !wrapper )
		{
			wrapper = document.createElement('div');
			wrapper.id = args.target;
			document.body.appendChild(wrapper);
		}
		//wrapper.innerHTML="";	// just to erase the initial HTML message
		wrapper.appendChild(rootWidget.domNode);

		return wrapper;
	},


	// Create a widget hierarchy from a JSON structure like
	// {widgetType: "LayoutContainer", params: { ... }, children: { ... } }
	// http://trac.dojotoolkit.org/browser/dijit/trunk/tests/layout/test_LayoutCode.html?rev=9408
	createWidgetHierarchy : function (widgetJson)
	{

		// setup input node
		node = document.createElement("div"),

		document.body.appendChild(node);	// necessary for tab contianer ???
		if(widgetJson.style){
			node.style.cssText = widgetJson.style;
		}
		if(widgetJson.innerHTML){
			node.innerHTML=widgetJson.innerHTML;
		}

		var widget;
		if ( widgetJson.widgetType == "dijit.form.Button" )
		{
			widget = new dijit.form.Button(widgetJson.params, node);
		}
		else if ( widgetJson.widgetType == "dijit.Toolbar" )
		{
			widget = new dijit.Toolbar(widgetJson.params, node);
		}
		else
		{
			var widgetType = widgetJson.widgetType.match(/^dijit.layout.(.+)$/)[1];
			widget = new dijit.layout[widgetType](widgetJson.params, node);
			widgetJson.widgetType = widgetType;
		}

		var widgetFramework = this;

		// add its children (recursively)
		if(widgetJson.children)
		{
			dojo.forEach( widgetJson.children,
				function(child)
				{

					widget.addChild(widgetFramework.createWidgetHierarchy(child));
				}
			);
		}


		widget.startup(); //TODO: this is required now, right?

		return widget;
	},

	// write out a menu of operations on each widget
	makeOperationTable : function ()
	{
		var html = "<table border=1>";
		dijit.registry.forEach(function(widget){
			html += "<tr><td>" + widget.id + "</td><td>";
			html += "<button onclick='removeFromParent(\"" + widget.id + "\");'> destroy </button> ";
			if(/Container/.test(widget.declaredClass)){
				html += "<button onclick='addChild(\"" + widget.id + "\");'> add a child </button> ";
			}
			html += "</td></tr>";
		});
		html += "</table>";

		if ( dojo.byId("operations") )
		{
			dojo.byId("operations").innerHTML = html;
		}
	},


   // remove a widget from it's parent and destroy it
   removeFromParent : function (widget)
   {
		widget = dijit.byId(widget);

		if(widget.parent)
		{
			widget.parent.removeChild(widget);
		}
		else if ( widget.parentNode )
		{
			widget.parentNode.removeChild(widget);
		}
		//else
		//{
		//}
		widget.destroy();

		// reset the operation table so this widget is no longer shown
		this.makeOperationTable();
	},

	// add a child to given widget
	addChild : function (widget)
	{
		widget = dijit.byId(widget);

		// setup input node
		var node = document.createElement("div");
		node.style.cssText = "height: 70px; width: 150px; overflow: auto; background: #cccccc; border: dotted black 2px;";	// necessary if parent is LayoutContainer
		// create the widget
		var alignments = ["top","bottom","left","right"];
		var hrefs = ["doc0.html", "doc1.html", "doc2.html"];
		var child = new dijit.layout.ContentPane(
			{
				title: "Widget " + this.cnt,	// necessary if parent is tab
				layoutAlign: alignments[this.cnt%4],	// necessary if parent is LayoutContainer
				executeScripts: true,
				href: hrefs[this.cnt%3]
			},
			node);
		this.cnt++;

		if(/AccordionContainer/.test(widget.declaredClass)){
			var pane = new dijit.layout.AccordionPane({
				title: "AccordionWidget " + this.cnt
			});
			pane.setContent(child);
			child = pane;
		}
		// add it to the parent
		widget.addChild(child);

		// reset the operation table so the new widget is shown
		makeOperationTable();
	},

	// show a widget
	show : function (widget){
		widget = dijit.byId(widget);
		widget.show();
	},

	// hide a widget
	hide : function (widget){
		widget = dijit.byId(widget);
		widget.hide();
	}

});  // end of WidgetFramework
