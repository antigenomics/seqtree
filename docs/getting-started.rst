Getting Started
===============

Install
-------

seqtree builds from source (C++17/20 core + pybind11). From a clone:

.. code-block:: fish

   bash setup.sh            # repo-local .venv + editable install
   bash setup.sh --tests    # also install pytest
   bash setup.sh --bench    # also install benchmark deps

Or directly with pip:

.. code-block:: console

   pip install -e .

First search
------------

.. code-block:: python

   import seqtree

   db = ["CASSLAPGATNEKLFF", "CASSLELGATNEKLFF", "CASSPQGATNEKLFF"]
   idx = seqtree.Index.build(db, alphabet="aa")   # "aa", "nt", or "iupac"

   p = seqtree.SearchParams(max_subs=2, engine="seqtm")
   for hit in idx.search("CASSLAPGATNEKLFF", p):
       print(hit.ref_id, hit.score, hit.n_subs, hit.n_ins, hit.n_dels)

Each :class:`~seqtree.Hit` is also tuple-unpackable:

.. code-block:: python

   ref_id, score, n_subs, n_ins, n_dels = hit

Batches in parallel
-------------------

``search_batch`` releases the GIL and runs a C++ thread pool over the queries:

.. code-block:: python

   results = idx.search_batch(queries, p, threads=0)   # 0 = all cores
   # results[i] is the hit list for queries[i]

Top hits, matrices, alignment
------------------------------

.. code-block:: python

   # k best hits
   top = idx.search_top("CASSLAPGATNEKLFF", p, k=5)

   # BLOSUM62-weighted budget (seqtrie)
   pm = seqtree.SearchParams(matrix="BLOSUM62", max_penalty=12, engine="seqtrie", gap_open=8)
   hits = idx.search("CASSLAPGATNEKLFF", pm)

   # alignment on demand (never computed during search)
   aln = idx.align(0, "CASSLELGATNEKLFF", p)
   print(aln.aligned_query, aln.aligned_ref, aln.ops)   # ops: M/S/I/D per column

Batch-vs-batch
--------------

For comparing two sets, ``pairwise_batch`` indexes the larger set automatically and
streams the smaller; results are always a-major:

.. code-block:: python

   pairs = seqtree.pairwise_batch(query_set, db_set, p, alphabet="aa")
   # pairs[i] are hits for query_set[i]; Hit.ref_id indexes db_set
