
from typing import Any, Dict, Iterable, List, Set, Tuple

import networkx as nx

class PathwayGraphFunctions:
    """Utility functions for pathway graph construction and analysis."""

    @staticmethod
    def make_pathway_graph(
        reactions_lists: Iterable[Tuple[str, str, str]],
        pathway_id: str,
    ) -> Tuple[
        nx.DiGraph,
        List[str],
        List[str],
        List[Dict[str, Any]],
        List[Dict[str, Any]],
    ]:
        """
        Build a directed pathway graph.

        Nodes are Reactome event stable IDs.
        Edges represent explicit BioPAX PathwayStep ordering.
        """

        graph = nx.DiGraph()

        for current, next_reaction, previous in reactions_lists:

            if not current:
                continue

            # Important: keep isolated reactions too.
            graph.add_node(current)

            if next_reaction:
                graph.add_edge(
                    current,
                    next_reaction,
                )

            if previous:
                graph.add_edge(
                    previous,
                    current,
                )

        starting_nodes = sorted(
            node
            for node in graph.nodes
            if graph.in_degree(node) == 0
        )

        ending_nodes = sorted(
            node
            for node in graph.nodes
            if graph.out_degree(node) == 0
        )

        edge_rows = [
            {
                "pathway_id": pathway_id,
                "source": source,
                "target": target,
            }
            for source, target in graph.edges
        ]

        node_rows = [
            {
                "pathway_id": pathway_id,
                "node": node,
                "is_start": node in starting_nodes,
                "is_end": node in ending_nodes,
                "in_degree": graph.in_degree(node),
                "out_degree": graph.out_degree(node),
                "in_cycle": False,
            }
            for node in graph.nodes
        ]

        cycle_nodes = {
            node
            for component in nx.strongly_connected_components(
                graph
            )
            if len(component) > 1
            for node in component
        }

        # Self-loops are cycles too.
        cycle_nodes.update(
            node
            for node in graph.nodes
            if graph.has_edge(node, node)
        )

        for row in node_rows:
            row["in_cycle"] = (
                row["node"] in cycle_nodes
            )

        return (
            graph,
            starting_nodes,
            ending_nodes,
            edge_rows,
            node_rows,
        )

    @staticmethod
    def remove_duplicates_order(lists: Iterable[List[Any]]) -> List[List[Any]]:
        """Remove duplicate lists while preserving order."""
        seen: Set[Tuple[Any, ...]] = set()
        result: List[List[Any]] = []
        for sub in lists:
            t = tuple(sub)
            if t not in seen:
                seen.add(t)
                result.append(list(sub))
        return result

    @staticmethod
    def is_subsequence(sub: List[Any], main: List[Any]) -> bool:
        """Return True if ``sub`` is a subsequence of ``main``."""
        it = iter(main)
        return all(x in it for x in sub)

    @staticmethod
    def find_non_subsequences(lists: Iterable[List[Any]]) -> List[List[Any]]:
        """Return those lists that are not subsequences of any other list."""
        lists = list(lists)
        non_subsequences: List[List[Any]] = []
        for i, sub in enumerate(lists):
            is_sub = False
            for j, main in enumerate(lists):
                if i != j and PathwayGraphFunctions.is_subsequence(sub, main):
                    is_sub = True
                    break
            if not is_sub:
                non_subsequences.append(sub)
        return non_subsequences

