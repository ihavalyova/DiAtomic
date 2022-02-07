{{ fullname | escape | underline}}

.. currentmodule:: {{ module }}

.. autoclass:: {{ objname }}
   :members:
   :undoc-members:
   :inherited-members:
   :show-inheritance:
   :private-members:

   {% block methods %}

   {% if methods %}
   .. rubric:: {{ ('Methods') }}

   .. autosummary::
   {% for item in members %}
      {%- if not item.startswith('__') %}
      ~{{ name }}.{{ item }}
      {%- endif -%}
   {%- endfor %}
   {% endif %}
   {% endblock %}

   {% block attributes %}
   {% if attributes %}
   .. rubric:: {{ ('Attributes') }}

   .. autosummary::
   {% for item in attributes %}
      ~{{ name }}.{{ item }}
   {%- endfor %}
   {% endif %}
   {% endblock %}

..
   :special-members:
