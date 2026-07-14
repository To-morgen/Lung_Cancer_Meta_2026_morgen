# Public Branch Workflow

This repository may be paired with local branches that contain unpublished
analysis notes or other private context. Removing those files in a later commit
does not remove them from Git history.

## History isolation rule

Never push a local analysis branch that has contained private handoffs, result
interpretation, machine paths, or manuscript-specific configuration.

Create a public feature branch from the current public base instead:

    git fetch origin
    git switch --create feature/<public-topic> origin/main

Then re-implement or cherry-pick only commits whose complete diffs have passed a
public-safety review. Do not merge the private branch.

## Local-only surfaces

Keep these paths untracked:

    configs_private/
    Notes/
    docs/HANDOFF_*.md
    data/
    results/
    figures/
    logs/

Private figure specifications belong under configs_private/figures/. Generic
schemas, synthetic fixtures, reusable builders, and tests may be public.

## Pre-push review

Review the full branch delta, not only the working tree:

    git status --short --branch
    git log --oneline origin/main..HEAD
    git diff --name-status origin/main...HEAD
    git diff --check origin/main...HEAD

Check tracked paths and content:

    git ls-files | grep -E '(^|/)(data|results|logs|configs_private|Notes|[.]snakemake)(/|$)' || true
    git ls-files 'docs/HANDOFF_*.md'
    git grep -n -E '/home/|/data[0-9]+/|token|password|api[_-]?key' origin/main..HEAD -- . || true

Push a reviewed feature branch first. Merge through the public review path
rather than pushing a private analysis branch directly to main.
